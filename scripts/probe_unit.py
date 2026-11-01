#!/usr/bin/env python3
"""
Probe exactly what aqi.in's GetDeviceLastHourData returns for one or more units.
Compares a suspect unit against a known-good one so we can see the difference.

    python3 probe_unit.py 2 6        # probe unit 2 and unit 6
    python3 probe_unit.py 2          # just unit 2 (defaults to comparing vs nothing)

Reads TOKEN / PARAMS / MONITORS from config.py.
"""
import sys, requests
from config import TOKEN, PARAMS, MONITORS

API = "https://airquality.aqi.in/api/v1/GetDeviceLastHourData"

def headers(serial, sid):
    return {"accept": "*/*", "authorization": f"bearer {TOKEN}",
            "origin": "https://dash.aqi.in", "referer": "https://dash.aqi.in/",
            "sensorid": str(sid), "serialno": str(serial),
            "user-agent": "Mozilla/5.0 (probe)"}

ids = sys.argv[1:] or ["2"]
for mid in ids:
    cfg = MONITORS.get(mid)
    if not cfg or not cfg.get("serialno"):
        print(f"\n[{mid}] no serial in config"); continue
    serial = cfg["serialno"]
    print(f"\n===== monitor {mid}   serial {serial}   role {cfg.get('role')} =====")
    for name, sid in PARAMS.items():
        if not sid:
            continue
        try:
            r = requests.get(API, headers=headers(serial, sid), timeout=15)
            j = r.json()
            d = j.get("data", {}) or {}
            vals = d.get("avrangearray") or []
            ts = d.get("timearray") or []
            if vals:
                print(f"  {name:>11} (sid {sid}): latest={vals[0]}  at {ts[0] if ts else '?'}  ({len(vals)} pts)")
            else:
                print(f"  {name:>11} (sid {sid}): EMPTY   data-keys={list(d.keys())}   http={r.status_code}")
        except Exception as e:
            print(f"  {name:>11} (sid {sid}): ERROR {e}")
