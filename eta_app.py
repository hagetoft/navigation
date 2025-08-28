#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import re
from datetime import datetime, timedelta
try:
    from zoneinfo import ZoneInfo  # Python 3.9+
except Exception:
    ZoneInfo = None

import streamlit as st
import pydeck as pdk

EARTH_RADIUS_NM = 3440.065  # nautical miles

def deg2rad(d): return d * math.pi / 180.0
def rad2deg(r): return r * 180.0 / math.pi

# ---- Parsing helpers ----
def parse_decimal(text: str) -> float:
    s = str(text).strip().replace(",", ".")
    return float(s)

def parse_dm(text: str, is_lat: bool) -> float:
    """Parse degrees + decimal minutes formats, e.g.:
    57°41.900'N, 57 41.900 N, N57 41.9, 11°59.000'E, -57 41.9
    Returns decimal degrees. is_lat picks valid range and hemisphere letters.
    """
    s = str(text).strip().upper()
    s = s.replace("º", "°").replace(",", ".")
    # Extract hemisphere letters if any
    hemi_letters = []
    for ch in "NSEW":
        if ch in s:
            hemi_letters.append(ch)
    sign = 1.0
    if hemi_letters:
        # Use the last letter if multiple were typed
        h = hemi_letters[-1]
        if is_lat and h in ("S",): sign = -1.0
        if is_lat and h in ("N",): sign = +1.0
        if (not is_lat) and h in ("W",): sign = -1.0
        if (not is_lat) and h in ("E",): sign = +1.0
    # Remove all non-digit separators except dot/plus/minus to find numbers
    s_numbers = re.sub(r"[NSEW]", " ", s)
    # Replace degree/minute/second markers with spaces
    s_numbers = re.sub(r"[°'\"]", " ", s_numbers)
    # Collapse spaces
    s_numbers = re.sub(r"\s+", " ", s_numbers).strip()
    nums = re.findall(r"[-+]?\d+(?:\.\d+)?", s_numbers)
    if len(nums) == 0:
        raise ValueError(f"Kunde inte tolka DM-koordinat: '{text}'")
    if len(nums) == 1:
        # Probably decimal degrees (fallback)
        val = float(nums[0])
        if hemi_letters and sign < 0:
            val = -abs(val)
        return val
    # Use first two numbers as deg + minutes
    deg = float(nums[0]); minutes = float(nums[1])
    if minutes < 0 or minutes >= 60:
        raise ValueError("Minuter måste vara mellan 0 och < 60")
    val = abs(deg) + minutes/60.0
    # If explicit sign in degrees part, respect it unless hemisphere given
    if deg < 0:
        sign = -1.0
    val *= sign
    # Validate ranges
    if is_lat and not (-90.0 <= val <= 90.0):
        raise ValueError("Latitud måste vara mellan -90 och 90")
    if (not is_lat) and not (-180.0 <= val <= 180.0):
        raise ValueError("Longitud måste vara mellan -180 och 180")
    return val

# ---- Core calculations ----
def distance_gc_nm(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    dphi = phi2 - phi1
    dlam = lam2 - lam1
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return EARTH_RADIUS_NM * c, c  # return central angle too

def initial_gc_bearing_deg(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    y = math.sin(lam2 - lam1) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(lam2 - lam1)
    theta = math.atan2(y, x)
    return (rad2deg(theta) + 360.0) % 360.0

def distance_rhumb_nm(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    dphi = phi2 - phi1
    dlam = abs(lam2 - lam1)
    if dlam > math.pi:
        dlam = 2*math.pi - dlam  # wrap across dateline the short way
    if phi1 != phi2:
        dpsi = math.log(math.tan(math.pi/4 + phi2/2) / math.tan(math.pi/4 + phi1/2))
    else:
        dpsi = 0.0
    q = dphi/dpsi if abs(dpsi) > 1e-12 else math.cos(phi1)
    dist = math.hypot(dphi, q * dlam) * EARTH_RADIUS_NM
    return dist

def rhumb_bearing_deg(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    dlam = lam2 - lam1
    if abs(dlam) > math.pi:
        dlam = dlam - math.copysign(2*math.pi, dlam)  # shortest wrap
    if phi1 != phi2:
        dpsi = math.log(math.tan(math.pi/4 + phi2/2) / math.tan(math.pi/4 + phi1/2))
    else:
        dpsi = 0.0
    if abs(dpsi) > 1e-12:
        theta = math.atan2(dlam, dpsi)
    else:
        theta = 0.0 if dlam >= 0 else math.pi
    return (rad2deg(theta) + 360.0) % 360.0

def travel_time(distance_nm, speed_knots):
    if speed_knots <= 0:
        raise ValueError("Speed in knots must be > 0.")
    hours = distance_nm / float(speed_knots)
    return timedelta(seconds=hours * 3600.0)

def fmt_td(td: timedelta):
    total_seconds = int(round(td.total_seconds()))
    h, r = divmod(total_seconds, 3600)
    m, s = divmod(r, 60)
    if h >= 24:
        d, h = divmod(h, 24)
        return f"{d}d {h:02d}:{m:02d}:{s:02d}"
    return f"{h:02d}:{m:02d}:{s:02d}"

# Great-circle interpolation for map path (slerp on the sphere)
def great_circle_points(lat1, lon1, lat2, lon2, n_points=64):
    phi1, lam1 = deg2rad(lat1), deg2rad(lon1)
    phi2, lam2 = deg2rad(lat2), deg2rad(lon2)
    x1 = math.cos(phi1) * math.cos(lam1)
    y1 = math.cos(phi1) * math.sin(lam1)
    z1 = math.sin(phi1)
    x2 = math.cos(phi2) * math.cos(lam2)
    y2 = math.cos(phi2) * math.sin(lam2)
    z2 = math.sin(phi2)
    dot = max(-1.0, min(1.0, x1*x2 + y1*y2 + z1*z2))
    omega = math.acos(dot)
    if omega == 0.0:
        return [[lon1, lat1], [lon2, lat2]]
    pts = []
    for i in range(n_points+1):
        f = i / n_points
        sin_1 = math.sin((1-f) * omega)
        sin_f = math.sin(f * omega)
        sin_o = math.sin(omega)
        x = (sin_1 * x1 + sin_f * x2) / sin_o
        y = (sin_1 * y1 + sin_f * y2) / sin_o
        z = (sin_1 * z1 + sin_f * z2) / sin_o
        lat = rad2deg(math.atan2(z, math.hypot(x, y)))
        lon = rad2deg(math.atan2(y, x))
        if lon > 180: lon -= 360
        if lon < -180: lon += 360
        pts.append([lon, lat])
    return pts

def rhumb_line_points(lat1, lon1, lat2, lon2, n_points=64):
    phi1, lam1 = deg2rad(lat1), deg2rad(lon1)
    phi2, lam2 = deg2rad(lat2), deg2rad(lon2)
    dlam = lam2 - lam1
    if abs(dlam) > math.pi:
        lam2 = lam2 - math.copysign(2*math.pi, dlam)
    def merc_y(phi): return math.log(math.tan(math.pi/4 + phi/2))
    y1 = merc_y(phi1)
    y2 = merc_y(phi2)
    pts = []
    for i in range(n_points+1):
        f = i / n_points
        lam = lam1 + f * (lam2 - lam1)
        y = y1 + f * (y2 - y1)
        phi = 2 * math.atan(math.exp(y)) - math.pi/2
        lat = rad2deg(phi)
        lon = rad2deg(lam)
        if lon > 180: lon -= 360
        if lon < -180: lon += 360
        pts.append([lon, lat])
    return pts

def zoom_by_distance_nm(distance_nm):
    if distance_nm < 10: return 11
    if distance_nm < 30: return 10
    if distance_nm < 100: return 8
    if distance_nm < 300: return 6
    if distance_nm < 1000: return 4
    return 3

def cumulative_distances_nm(points):
    cum = [0.0]
    total = 0.0
    for i in range(1, len(points)):
        lon1, lat1 = points[i-1]
        lon2, lat2 = points[i]
        seg_nm, _ = distance_gc_nm(lat1, lon1, lat2, lon2)
        total += seg_nm
        cum.append(total)
    return cum

def generate_gpx(points, start, dest, method_name, speed_knots, dt_dep=None, include_times=False, route_name="ETA Route"):
    now_iso = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    rtepts = [f'<rtept lat="{lat:.6f}" lon="{lon:.6f}"></rtept>' for lon, lat in points]
    trkpts = []
    if include_times and dt_dep is not None and speed_knots > 0:
        cum = cumulative_distances_nm(points)
        for i, (lon, lat) in enumerate(points):
            hours = cum[i] / float(speed_knots)
            t = dt_dep + timedelta(hours=hours)
            t_iso = t.astimezone(ZoneInfo("UTC")).strftime("%Y-%m-%dT%H:%M:%SZ") if ZoneInfo else t.strftime("%Y-%m-%dT%H:%M:%SZ")
            trkpts.append(f'<trkpt lat="{lat:.6f}" lon="{lon:.6f}"><time>{t_iso}</time></trkpt>')
    else:
        for lon, lat in points:
            trkpts.append(f'<trkpt lat="{lat:.6f}" lon="{lon:.6f}"></trkpt>')
    wpt_start = f'<wpt lat="{start[0]:.6f}" lon="{start[1]:.6f}"><name>Start</name></wpt>'
    wpt_dest  = f'<wpt lat="{dest[0]:.6f}" lon="{dest[1]:.6f}"><name>Destination</name></wpt>'
    gpx = f'''<?xml version="1.0" encoding="UTF-8"?>
<gpx version="1.1" creator="ETA Streamlit" xmlns="http://www.topografix.com/GPX/1/1"
     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     xsi:schemaLocation="http://www.topografix.com/GPX/1/1 http://www.topografix.com/GPX/1/1/gpx.xsd">
  <metadata>
    <name>{route_name}</name>
    <time>{now_iso}</time>
    <desc>Method: {method_name}; Speed: {speed_knots:.2f} kn</desc>
  </metadata>
  {wpt_start}
  {wpt_dest}
  <rte>
    <name>{route_name}</name>
    <desc>{method_name} route</desc>
    {''.join(rtepts)}
  </rte>
  <trk>
    <name>{route_name}</name>
    <type>{method_name}</type>
    <trkseg>
      {''.join(trkpts)}
    </trkseg>
  </trk>
</gpx>'''
    return gpx

# ---- UI ----
st.set_page_config(page_title="ETA-kalkylator", layout="wide")
st.title("ETA-kalkylator med karta, GPX och DM-inmatning")

with st.form("inputs"):
    fmt = st.radio("Koordinatformat", ["Decimalgrader", "Grader + minuter (DM.m)"])
    c1, c2 = st.columns(2)
    if fmt == "Decimalgrader":
        with c1:
            st.subheader("Start")
            lat1_val = st.number_input("Lat (decimal)", value=57.6900, step=0.0001, format="%.6f")
            lon1_val = st.number_input("Lon (decimal)", value=11.9833, step=0.0001, format="%.6f")
        with c2:
            st.subheader("Mål")
            lat2_val = st.number_input("Lat (decimal)", value=57.7000, step=0.0001, format="%.6f", key="lat2")
            lon2_val = st.number_input("Lon (decimal)", value=11.6000, step=0.0001, format="%.6f", key="lon2")
    else:
        with c1:
            st.subheader("Start (DM.m)")
            lat1_text = st.text_input("Lat (t.ex. 57°41.900'N eller 57 41.900 N)", value="57°41.400'N")
            lon1_text = st.text_input("Lon (t.ex. 11°59.000'E eller 11 59.000 E)", value="11°59.000'E")
        with c2:
            st.subheader("Mål (DM.m)")
            lat2_text = st.text_input("Lat", value="57°42.000'N", key="lat2t")
            lon2_text = st.text_input("Lon", value="11°36.000'E", key="lon2t")
    c3, c4, c5 = st.columns(3)
    with c3:
        speed = st.number_input("Medelhastighet (knop)", min_value=0.1, value=20.0, step=0.1, format="%.1f")
    with c4:
        method = st.selectbox("Metod", options=["Great-circle", "Rhumb line"], index=0)
    with c5:
        use_depart = st.checkbox("Ange avresetid", value=False)
    include_times = False
    if use_depart:
        dcol, tcol = st.columns(2)
        with dcol:
            dep_date = st.date_input("Datum", value=datetime.now(ZoneInfo("Europe/Stockholm")).date() if ZoneInfo else datetime.now().date())
        with tcol:
            dep_time = st.time_input("Tid", value=datetime.now(ZoneInfo("Europe/Stockholm")).time() if ZoneInfo else datetime.now().time())
        include_times = st.checkbox("Inkludera tidsstämplar i GPX", value=True)
    route_name = st.text_input("Rutt-namn (för GPX)", value="ETA Route")
    submit = st.form_submit_button("Beräkna")

if submit:
    try:
        # Parse coords
        if fmt == "Decimalgrader":
            lat1, lon1, lat2, lon2 = float(lat1_val), float(lon1_val), float(lat2_val), float(lon2_val)
        else:
            lat1 = parse_dm(lat1_text, is_lat=True)
            lon1 = parse_dm(lon1_text, is_lat=False)
            lat2 = parse_dm(lat2_text, is_lat=True)
            lon2 = parse_dm(lon2_text, is_lat=False)

        if method == "Great-circle":
            dist_nm, _ = distance_gc_nm(lat1, lon1, lat2, lon2)
            bearing = initial_gc_bearing_deg(lat1, lon1, lat2, lon2)
            path = great_circle_points(lat1, lon1, lat2, lon2, n_points=256)
            method_name = "Great-circle"
        else:
            dist_nm = distance_rhumb_nm(lat1, lon1, lat2, lon2)
            bearing = rhumb_bearing_deg(lat1, lon1, lat2, lon2)
            path = rhumb_line_points(lat1, lon1, lat2, lon2, n_points=256)
            method_name = "Rhumb line"

        tt = travel_time(dist_nm, speed)
        if use_depart:
            tz = ZoneInfo("Europe/Stockholm") if ZoneInfo else None
            dt_dep = datetime.combine(dep_date, dep_time, tzinfo=tz)
        else:
            dt_dep = datetime.now(ZoneInfo("Europe/Stockholm")) if ZoneInfo else datetime.now()
        dt_eta = dt_dep + tt

        # Results
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Avstånd", f"{dist_nm:.2f} NM")
        m2.metric("Initial kurs", f"{bearing:.1f}°")
        m3.metric("Restid", fmt_td(tt))
        tzstr = (dt_eta.tzname() or "") if hasattr(dt_eta, "tzname") else ""
        m4.metric("ETA", dt_eta.strftime("%Y-%m-%d %H:%M") + (f" {tzstr}" if tzstr else ""))

        # Map
        avg_lat = (lat1 + lat2) / 2.0
        avg_lon = (lon1 + lon2) / 2.0
        zoom = zoom_by_distance_nm(dist_nm)

        path_layer = pdk.Layer(
            "PathLayer",
            data=[{"path": path}],
            get_path="path",
            get_width=3,
            width_min_pixels=2,
            get_color=[30, 144, 255, 200],
            pickable=False,
        )
        points = [
            {"name": "Start", "lon": lon1, "lat": lat1},
            {"name": "Mål", "lon": lon2, "lat": lat2},
        ]
        point_layer = pdk.Layer(
            "ScatterplotLayer",
            data=points,
            get_position="[lon, lat]",
            get_radius=100,
            radius_min_pixels=6,
            get_fill_color=[255, 99, 71, 200],
            pickable=True,
        )
        tooltip = {"text": "{name}\n{lat}, {lon}"}
        view_state = pdk.ViewState(latitude=avg_lat, longitude=avg_lon, zoom=zoom)
        st.pydeck_chart(pdk.Deck(layers=[path_layer, point_layer], initial_view_state=view_state, tooltip=tooltip, map_style=None))

        # GPX export
        gpx_str = generate_gpx(
            points=path,
            start=(lat1, lon1),
            dest=(lat2, lon2),
            method_name=method_name,
            speed_knots=speed,
            dt_dep=dt_dep,
            include_times=include_times,
            route_name=route_name,
        )
        st.download_button(
            "⬇️ Ladda ner GPX",
            data=gpx_str.encode("utf-8"),
            file_name=f"{route_name.replace(' ', '_')}.gpx",
            mime="application/gpx+xml",
        )

        with st.expander("Exempel på DM-inmatning"):
            st.markdown("""
            - `57°41.900'N`  (lat), `11°59.000'E` (lon)  
            - `57 41.900 N`, `11 59.000 E`  
            - `N57 41.900`, `E011 59.000`  
            - Negativa tecken funkar också: `-57 41.900` (S), `-011 59.000` (W)
            """)
    except Exception as e:
        st.error(f"Fel: {e}")

st.caption("Obs: geometrisk ETA utan vind/ström/hinder. Nästa steg kan vara vektorström, waypoints och sjökortsdata (S57/S101). ")
