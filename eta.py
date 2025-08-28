#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ETA calculator (ASCII-safe) between two coordinates with speed in knots.

Examples:
  python eta.py --start "57.6900,11.9833" --dest "57.7000,11.6000" --speed 25
  python eta.py --start 57.69 11.9833 --dest 57.70 11.60 --speed 20 --method rhumb --depart "2025-08-28T14:30"
"""
import math
import argparse
from datetime import datetime, timedelta
try:
    from zoneinfo import ZoneInfo  # Python 3.9+
except Exception:
    ZoneInfo = None

EARTH_RADIUS_NM = 3440.065  # nautical miles

def deg2rad(d): return d * math.pi / 180.0
def rad2deg(r): return r * 180.0 / math.pi

def parse_latlon(a, b=None):
    if b is None:
        if isinstance(a, str) and ',' in a:
            lat_s, lon_s = a.split(',', 1)
            return float(lat_s.strip()), float(lon_s.strip())
        raise ValueError("If --start/--dest is one value, it must be 'lat,lon'")
    return float(a), float(b)

# Great-circle (haversine)
def distance_gc_nm(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    dphi = phi2 - phi1
    dlam = lam2 - lam1
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return EARTH_RADIUS_NM * c

def initial_gc_bearing_deg(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    y = math.sin(lam2 - lam1) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(lam2 - lam1)
    theta = math.atan2(y, x)
    return (rad2deg(theta) + 360.0) % 360.0

# Rhumb line (loxodrome)
def distance_rhumb_nm(lat1, lon1, lat2, lon2):
    phi1, lam1, phi2, lam2 = map(deg2rad, [lat1, lon1, lat2, lon2])
    dphi = phi2 - phi1
    dlam = abs(lam2 - lam1)
    if dlam > math.pi:
        dlam = 2*math.pi - dlam  # dateline wrap shortest way
    # Mercator latitude difference
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

def parse_departure(depart_str):
    if not depart_str:
        return datetime.now(ZoneInfo("Europe/Stockholm")) if ZoneInfo else datetime.now()
    # ISO 8601 (optionally with TZ)
    try:
        dt = datetime.fromisoformat(depart_str)
        if dt.tzinfo is None and ZoneInfo:
            dt = dt.replace(tzinfo=ZoneInfo("Europe/Stockholm"))
        return dt
    except ValueError:
        pass
    # Fallback: YYYY-MM-DD HH:MM
    try:
        dt = datetime.strptime(depart_str, "%Y-%m-%d %H:%M")
        if ZoneInfo:
            dt = dt.replace(tzinfo=ZoneInfo("Europe/Stockholm"))
        return dt
    except ValueError:
        raise ValueError("Could not parse --depart. Use ISO 8601 like '2025-08-28T14:30'")

def main():
    p = argparse.ArgumentParser(description="ETA between two coordinates at a given average speed in knots.")
    p.add_argument("--start", nargs="+", required=True, help="Start as 'lat,lon' or two numbers 'lat lon'")
    p.add_argument("--dest", nargs="+", required=True, help="Destination as 'lat,lon' or two numbers 'lat lon'")
    p.add_argument("--speed", type=float, required=True, help="Average speed in knots (NM/h)")
    p.add_argument("--depart", type=str, default=None, help="Departure time (ISO 8601). Default: now (Europe/Stockholm)")
    p.add_argument("--method", choices=["gc", "rhumb"], default="gc", help="Distance method: gc=great-circle, rhumb=rhumb line")
    args = p.parse_args()

    # Parse coords
    lat1, lon1 = parse_latlon(*args.start) if len(args.start) == 1 else parse_latlon(args.start[0], args.start[1])
    lat2, lon2 = parse_latlon(*args.dest)  if len(args.dest)  == 1 else parse_latlon(args.dest[0],  args.dest[1])

    # Compute
    if args.method == "gc":
        dist_nm = distance_gc_nm(lat1, lon1, lat2, lon2)
        brg_deg = initial_gc_bearing_deg(lat1, lon1, lat2, lon2)
        method_name = "Great-circle"
    else:
        dist_nm = distance_rhumb_nm(lat1, lon1, lat2, lon2)
        brg_deg = rhumb_bearing_deg(lat1, lon1, lat2, lon2)
        method_name = "Rhumb line"

    tt = travel_time(dist_nm, args.speed)
    depart_dt = parse_departure(args.depart)
    eta_dt = depart_dt + tt

    # Output
    def fmt_td(td: timedelta):
        total_seconds = int(round(td.total_seconds()))
        h, r = divmod(total_seconds, 3600)
        m, s = divmod(r, 60)
        if h >= 24:
            d, h = divmod(h, 24)
            return f"{d}d {h:02d}:{m:02d}:{s:02d}"
        return f"{h:02d}:{m:02d}:{s:02d}"

    print(f"Method: {method_name}")
    print(f"Start (lat, lon): {lat1:.6f}, {lon1:.6f}")
    print(f"Dest  (lat, lon): {lat2:.6f}, {lon2:.6f}")
    print(f"Initial course:   {brg_deg:0.1f}°")
    print(f"Distance:         {dist_nm:.2f} NM")
    print(f"Speed:            {args.speed:.2f} kn")
    print(f"Travel time:      {fmt_td(tt)}")
    print("Departure:        " + depart_dt.strftime('%Y-%m-%d %H:%M:%S %Z'))
    print("ETA:               " + eta_dt.strftime('%Y-%m-%d %H:%M:%S %Z'))

if __name__ == "__main__":
    main()
