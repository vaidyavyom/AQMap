#!/usr/bin/env python3
"""
Fetch ERA5 hourly 10m wind (u, v components) + 850 hPa wind for Chandrapur.

Setup (one-time):
  1. Register at https://cds.climate.copernicus.eu
  2. Go to your profile page and copy your API key
  3. Create ~/.cdsapirc with:
       url: https://cds.climate.copernicus.eu/api
       key: YOUR-API-KEY-HERE

Install:
  pip3 install cdsapi --break-system-packages

Run:
  python3 fetch_era5_wind.py

Outputs:
  era5_wind_chandrapur.csv  — hourly wind speed + direction at 10m and 850 hPa
                              covering Aug 11 – Oct 6 2026 (full deployment window)
"""

import cdsapi, math, csv
from datetime import datetime, timedelta, timezone
from pathlib import Path

# ── bounding box around Chandrapur network (0.25° ERA5 grid) ──────────────
LAT = 19.97    # roughly centre of monitor cluster
LON = 79.29
# ERA5 area: N/W/S/E — one grid cell buffer each side
AREA = [20.25, 79.00, 19.75, 79.50]   # N, W, S, E

# ── deployment window ──────────────────────────────────────────────────────
START = datetime(2026, 8, 11)
END   = datetime(2026, 10, 6, 23, 59)

OUT = "era5_wind_chandrapur.csv"


def uv_to_dir_speed(u, v):
    """Meteorological wind direction (where wind comes FROM) and speed."""
    speed = math.hypot(u, v)
    # math.atan2 gives direction wind blows TO; add 180 for FROM
    direction = (270 - math.degrees(math.atan2(v, u))) % 360
    return round(direction, 1), round(speed, 2)


def build_month_list(start, end):
    """List of (year, month) tuples spanning the date range."""
    months = []
    y, m = start.year, start.month
    while (y, m) <= (end.year, end.month):
        months.append((y, m))
        m += 1
        if m > 12:
            m = 1
            y += 1
    return months


def fetch_era5(pressure_level=False):
    c = cdsapi.Client()
    months = build_month_list(START, END)

    all_files = []
    for year, month in months:
        # days in this month that fall within our window
        import calendar
        _, ndays = calendar.monthrange(year, month)
        days = [str(d) for d in range(1, ndays + 1)]

        fname = f"era5_raw_{year}_{month:02d}{'_pl' if pressure_level else '_sfc'}.nc"
        if Path(fname).exists():
            print(f"  {fname} already cached, skipping download")
            all_files.append(fname)
            continue

        print(f"  Downloading {year}-{month:02d} ({'pressure level' if pressure_level else 'surface'})...")

        if pressure_level:
            c.retrieve(
                "reanalysis-era5-pressure-levels",
                {
                    "product_type": "reanalysis",
                    "variable": ["u_component_of_wind", "v_component_of_wind"],
                    "pressure_level": ["850"],
                    "year": str(year),
                    "month": f"{month:02d}",
                    "day": days,
                    "time": [f"{h:02d}:00" for h in range(24)],
                    "area": AREA,
                    "format": "netcdf",
                },
                fname,
            )
        else:
            c.retrieve(
                "reanalysis-era5-single-levels",
                {
                    "product_type": "reanalysis",
                    "variable": [
                        "10m_u_component_of_wind",
                        "10m_v_component_of_wind",
                    ],
                    "year": str(year),
                    "month": f"{month:02d}",
                    "day": days,
                    "time": [f"{h:02d}:00" for h in range(24)],
                    "area": AREA,
                    "format": "netcdf",
                },
                fname,
            )
        all_files.append(fname)

    return all_files


def extract_to_csv(sfc_files, pl_files):
    """Parse downloaded NetCDF files and write a clean CSV."""
    import netCDF4 as nc
    import numpy as np

    rows = {}  # datetime_utc -> {u10, v10, u850, v850}

    def nearest_idx(arr, val):
        return int(np.argmin(np.abs(arr - val)))

    for fpath in sfc_files:
        ds = nc.Dataset(fpath)
        tvar = "valid_time" if "valid_time" in ds.variables else "time"
        times = nc.num2date(ds.variables[tvar][:], ds.variables[tvar].units)
        latvar = "latitude" if "latitude" in ds.variables else "lat"
        lonvar = "longitude" if "longitude" in ds.variables else "lon"
        lats  = ds.variables[latvar][:]
        lons  = ds.variables[lonvar][:]
        li = nearest_idx(lats, LAT)
        lo = nearest_idx(lons, LON)

        u_var = "u10" if "u10" in ds.variables else "u_component_of_wind"
        v_var = "v10" if "v10" in ds.variables else "v_component_of_wind"

        for i, t in enumerate(times):
            dt = datetime(t.year, t.month, t.day, t.hour, tzinfo=timezone.utc)
            if not (START.replace(tzinfo=timezone.utc) <= dt <= END.replace(tzinfo=timezone.utc)):
                continue
            rows.setdefault(dt, {})
            rows[dt]["u10"] = float(ds.variables[u_var][i, li, lo])
            rows[dt]["v10"] = float(ds.variables[v_var][i, li, lo])
        ds.close()

    for fpath in pl_files:
        ds = nc.Dataset(fpath)
        tvar2 = "valid_time" if "valid_time" in ds.variables else "time"
        times = nc.num2date(ds.variables[tvar2][:], ds.variables[tvar2].units)
        latvar2 = "latitude" if "latitude" in ds.variables else "lat"
        lonvar2 = "longitude" if "longitude" in ds.variables else "lon"
        lats  = ds.variables[latvar2][:]
        lons  = ds.variables[lonvar2][:]
        li = nearest_idx(lats, LAT)
        lo = nearest_idx(lons, LON)
        # pressure axis
        plev  = ds.variables["pressure_level"][:] if "pressure_level" in ds.variables else ds.variables["level"][:]
        pi = nearest_idx(plev, 850)

        u_var = "u" if "u" in ds.variables else "u_component_of_wind"
        v_var = "v" if "v" in ds.variables else "v_component_of_wind"

        for i, t in enumerate(times):
            dt = datetime(t.year, t.month, t.day, t.hour, tzinfo=timezone.utc)
            if not (START.replace(tzinfo=timezone.utc) <= dt <= END.replace(tzinfo=timezone.utc)):
                continue
            rows.setdefault(dt, {})
            rows[dt]["u850"] = float(ds.variables[u_var][i, pi, li, lo])
            rows[dt]["v850"] = float(ds.variables[v_var][i, pi, li, lo])
        ds.close()

    # write CSV
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "datetime_utc", "datetime_ist",
            "u10_ms", "v10_ms", "wind_dir_10m", "wind_speed_10m_ms",
            "u850_ms", "v850_ms", "wind_dir_850hpa", "wind_speed_850hpa_ms",
        ])
        for dt in sorted(rows):
            r = rows[dt]
            ist = dt + timedelta(hours=5, minutes=30)
            d10, s10 = uv_to_dir_speed(r.get("u10", 0), r.get("v10", 0))
            d85, s85 = uv_to_dir_speed(r.get("u850", 0), r.get("v850", 0))
            w.writerow([
                dt.isoformat(), ist.replace(tzinfo=None).isoformat(),
                round(r.get("u10", ""), 3), round(r.get("v10", ""), 3), d10, s10,
                round(r.get("u850", ""), 3), round(r.get("v850", ""), 3), d85, s85,
            ])

    print(f"\nWrote {len(rows)} hourly rows → {OUT}")


if __name__ == "__main__":
    print("Fetching ERA5 surface wind (10m)...")
    sfc = fetch_era5(pressure_level=False)

    print("Fetching ERA5 pressure-level wind (850 hPa, ~1500m — better for tall stacks)...")
    pl  = fetch_era5(pressure_level=True)

    print("Extracting to CSV...")
    try:
        extract_to_csv(sfc, pl)
    except ImportError:
        print("\nnetCDF4 not installed. Run:")
        print("  pip3 install netCDF4 --break-system-packages")
        print("Then re-run this script.")
