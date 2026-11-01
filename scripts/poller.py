#!/usr/bin/env python3
"""
Chandrapur AQ network poller (aqi.in / GetDeviceLastHourData).
Reads settings from config.py, pulls each monitor's latest reading for every
parameter in PARAMS, applies calibration to PM2.5, writes latest.json and
appends history.csv.

  pip3 install requests
  python3 poller.py

Your token, monitors, and parameters live in config.py — edit that, not this.
"""
import json, time, csv, os, datetime, traceback, math, subprocess
import requests
from config import (TOKEN, PARAMS, INTERVAL_SEC, MONITORS, CALIB,
                    RH_CORRECTION, RH_KAPPA, RH_CAP)
try:
    from config import OPENAQ_KEY, OPENAQ_REFS
except ImportError:
    OPENAQ_KEY, OPENAQ_REFS = "", []

OUT_JSON, OUT_HISTORY, HEARTBEAT = "latest.json", "history.csv", "poller.log"
API_URL = "https://airquality.aqi.in/api/v1/GetDeviceLastHourData"
OPENAQ_URL = "https://api.openaq.org/v3"
DEBUG_DUMP = False   # set True to log raw responses while troubleshooting


def headers_for(serialno, sensorid):
    return {
        "accept": "*/*",
        "authorization": f"bearer {TOKEN}",
        "origin": "https://dash.aqi.in",
        "referer": "https://dash.aqi.in/",
        "sensorid": str(sensorid),
        "serialno": str(serialno),
        "user-agent": "Mozilla/5.0 (poller)",
    }

def fetch_value(serialno, sensorid):
    """latest value of one parameter for one device -> (value, ts) or (None, None)"""
    r = requests.get(API_URL, headers=headers_for(serialno, sensorid), timeout=15)
    r.raise_for_status()
    j = r.json()
    if DEBUG_DUMP:
        log("RAW %s/%s: %s" % (serialno, sensorid, json.dumps(j)[:300]))
    d = j.get("data", {}) or {}
    vals = d.get("avrangearray") or []
    ts = d.get("timearray") or []
    if not vals:
        return None, None
    return float(vals[0]), (ts[0] if ts else None)

def fetch_device(serialno):
    """all configured parameters for one device"""
    out = {}
    ts_seen = None
    for name, sid in PARAMS.items():
        if not sid:
            continue
        v, ts = fetch_value(serialno, sid)
        if v is not None:
            out[name] = v
            ts_seen = ts_seen or ts
    out["ts"] = ts_seen or datetime.datetime.now().isoformat()
    return out

def apply_calibration(mid, raw, rh=None):
    """monitor-to-monitor harmonization (+ optional reference-fit linear RH term)."""
    c = CALIB.get(mid, {})
    val = c.get("b0", 0.0) + c.get("b_raw", 1.0) * raw
    b_rh = c.get("b_rh", 0.0)          # 0 until a reference-fit RH term exists
    if b_rh and rh is not None:
        val += b_rh * rh
    return val

def rh_correct(pm, rh):
    """interim kappa-Kohler hygroscopic-growth correction (RH in %)."""
    if not RH_CORRECTION or RH_KAPPA <= 0 or rh is None:
        return pm
    aw = min(float(rh), RH_CAP) / 100.0
    if aw >= 0.999:
        aw = RH_CAP / 100.0
    gf = 1.0 + RH_KAPPA * aw / (1.0 - aw)   # volume growth factor
    return pm / gf

# ---- OpenAQ v3 reference stations (live coverage nodes) ----
_OAQ_SMAP = {}   # location_id -> {sensorsId: parameter_name}, cached
_OAQ_PARAM = {"pm25": "pm25", "pm10": "pm10", "relativehumidity": "humidity",
              "temperature": "temperature"}

def _oaq_get(path):
    r = requests.get(f"{OPENAQ_URL}{path}", headers={"X-API-Key": OPENAQ_KEY}, timeout=20)
    r.raise_for_status()
    return r.json()

def openaq_sensor_map(loc):
    res = _oaq_get(f"/locations/{loc}")["results"][0]["sensors"]
    return {s["id"]: s["parameter"]["name"] for s in res}

def fetch_openaq_ref(ref):
    """latest reference-grade reading for one OpenAQ location -> dict or None."""
    loc = ref["location_id"]
    if loc not in _OAQ_SMAP:
        _OAQ_SMAP[loc] = openaq_sensor_map(loc)
    smap = _OAQ_SMAP[loc]
    d = {}
    for row in _oaq_get(f"/locations/{loc}/latest").get("results", []):
        p = _OAQ_PARAM.get(smap.get(row.get("sensorsId")))
        if p and row.get("value") is not None:
            d[p] = float(row["value"])
    if "pm25" not in d:
        return None
    d["ts"] = datetime.datetime.now().isoformat()
    d["source"] = "openaq-ref"
    return d

def log(msg):
    line = f"{datetime.datetime.now().isoformat()}  {msg}"
    print(line)
    with open(HEARTBEAT, "a") as f:
        f.write(line + "\n")

# ---- coordinates + surroundings estimate (for units with a dead temp/RH board) ----
COORDS = {
    "1": (19.98062, 79.23735), "2": (19.97237, 79.28563), "3": (19.98250, 79.30010),
    "4": (19.96839, 79.31393), "5": (19.96024, 79.31746), "6": (19.97239, 79.30880),
    "7": (20.00482, 79.30233), "8": (19.94894, 79.29743), "9": (19.93473, 79.32349),
    "10": (19.95375, 79.32983), "R1": (19.9629, 79.298714),
}
def kmdist(la1, lo1, la2, lo2):
    dy = (la2 - la1) * 111.32
    dx = (lo2 - lo1) * 111.32 * math.cos(math.radians(la1))
    return math.hypot(dx, dy)

def idw_estimate(target_id, field, monitors):
    """estimate a met field at target_id from the other monitors that have it (IDW, p2)."""
    if target_id not in COORDS:
        return None
    tla, tlo = COORDS[target_id]
    num = den = 0.0
    for mid, d in monitors.items():
        if mid == target_id or mid not in COORDS:
            continue
        v = d.get(field)
        if v is None:
            continue
        dist = max(kmdist(tla, tlo, *COORDS[mid]), 0.05)
        w = 1.0 / dist**2
        num += w * v; den += w
    return num / den if den else None

def poll_once():
    out = {"generated": datetime.datetime.now().isoformat(), "monitors": {}}
    rows = []

    # pass 1: fetch every active device
    fetched = {}
    polled = []
    for mid, cfg in MONITORS.items():
        if cfg.get("role") == "pending" or not cfg.get("serialno"):
            continue
        polled.append(mid)
        try:
            d = fetch_device(cfg["serialno"])
        except Exception as e:
            log(f"WARN monitor {mid} fetch failed: {e}")
            continue
        if "pm25" in d:
            fetched[mid] = d

    # network-median RH (skip missing and stuck-at-100 sensors, e.g. unit 8)
    rhs = sorted(d["humidity"] for d in fetched.values()
                 if d.get("humidity") is not None and d["humidity"] < 99)
    med_rh = rhs[len(rhs)//2] if rhs else None

    # pass 2: calibrate (monitor-to-monitor) then apply interim RH correction
    for mid, d in fetched.items():
        raw_pm = d["pm25"]
        rh = d.get("humidity")
        rh_used = med_rh if (rh is None or rh >= 99) else rh    # dead RH -> network median
        harm = apply_calibration(mid, raw_pm)                   # monitor-to-monitor
        d["pm25"] = round(rh_correct(harm, rh_used), 1)         # + interim RH correction
        d["pm25_raw"] = raw_pm
        # a faulty env board (e.g. unit 8: temp 146C, RH stuck 100) must not poison the
        # map. if temperature is impossible, distrust humidity from the same board too;
        # both get filled from the surroundings below.
        temp = d.get("temperature")
        if temp is not None and not (-10 <= temp <= 55):
            d["temperature"] = None
            d["humidity"] = None
        elif d.get("humidity") is not None and not (0 < d["humidity"] <= 100):
            d["humidity"] = None
        out["monitors"][mid] = d

    # OpenAQ reference nodes (reference-grade; no calibration applied)
    if OPENAQ_KEY:
        for ref in OPENAQ_REFS:
            try:
                d = fetch_openaq_ref(ref)
            except Exception as e:
                log(f"WARN openaq ref {ref['id']} fetch failed: {e}")
                continue
            if not d:
                continue
            out["monitors"][ref["id"]] = d

    # fill dead/missing temp & humidity from the surrounding network (IDW), e.g. unit 8
    for mid, d in out["monitors"].items():
        for field in ("temperature", "humidity"):
            if d.get(field) is None:
                est = idw_estimate(mid, field, out["monitors"])
                if est is not None:
                    d[field] = round(est, 1)
                    d.setdefault("estimated", []).append(field)

    # build history rows from the finalized monitors (raw vs calibrated pm25)
    for mid, d in out["monitors"].items():
        rows.append([out["generated"], mid, d.get("pm25_raw", d["pm25"]), d["pm25"],
                     d.get("pm10", ""), d.get("temperature", ""), d.get("humidity", ""),
                     d.get("ts", "")])

    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(out, f)
    os.replace(tmp, OUT_JSON)

    new = not os.path.exists(OUT_HISTORY)
    with open(OUT_HISTORY, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["poll_time","monitor","raw_pm25","calibrated_pm25","pm10","temp","humidity","device_ts"])
        w.writerows(rows)
    nodata = [m for m in polled if m not in fetched]
    log(f"OK wrote {len(out['monitors'])} monitors"
        + (f"  | polled-but-no-data: {','.join(nodata)}" if nodata else ""))

GIT_PUSH_EVERY = 5        # push to GitHub every N polls (≈5 min at 60s interval)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def git_push():
    """Push latest.json + history.csv to GitHub Pages."""
    try:
        subprocess.run(
            ["git", "-C", ROOT, "add", "latest.json", "history.csv"],
            check=True, capture_output=True)
        ts = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
        result = subprocess.run(
            ["git", "-C", ROOT, "commit", "-m", f"data: {ts}"],
            capture_output=True, text=True)
        if "nothing to commit" in result.stdout:
            return   # no new data, skip push
        subprocess.run(
            ["git", "-C", ROOT, "push"],
            check=True, capture_output=True)
        log("pushed latest.json + history.csv to GitHub")
    except subprocess.CalledProcessError as e:
        log(f"WARN git push failed: {e.stderr.decode().strip()[:200]}")
    except Exception as e:
        log(f"WARN git push error: {e}")

def main():
    log(f"poller starting, interval={INTERVAL_SEC}s, params={[k for k,v in PARAMS.items() if v]}")
    poll_count = 0
    while True:
        try:
            poll_once()
            poll_count += 1
            if poll_count % GIT_PUSH_EVERY == 0:
                git_push()
        except Exception:
            log("ERROR:\n" + traceback.format_exc())
        time.sleep(INTERVAL_SEC)

if __name__ == "__main__":
    main()
