#!/usr/bin/env python3
"""
fig7_dispersion.py  —  Observed Δ PM2.5 vs Gaussian plume model
Compares per-station wind-sector PM2.5 differences against a Gaussian
ground-level concentration curve calibrated to CSTPS stack parameters
(CSTPS_dispersion_brief.pdf: Hs=275m, Vs=22m/s, d=6.5m, Ts=400K).

Formula: full 2D centerline C ∝ 1/(σy·σz)·exp(-He²/2σz²)
He is fit as a free parameter; Briggs estimate shown for reference.

Run from AQMap root:  python3 scripts/fig7_dispersion.py
Output:               figures/fig7_dispersion.png
"""

import numpy as np
import pandas as pd
from scipy import stats, optimize
import math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT       = os.path.dirname(SCRIPT_DIR)
DATA_DIR   = os.path.join(ROOT, "data", "sensors")
ERA5_PATH  = os.path.join(ROOT, "data", "era5", "era5_wind_chandrapur.csv")
FIG_OUT    = os.path.join(ROOT, "figures", "fig7_dispersion.png")

# ── Station geometry ───────────────────────────────────────────────────────
PLANT = (20.00837, 79.29206)
COORDS = {
    "1": (19.98062,79.23735),"2": (19.97237,79.28563),
    "3": (19.98250,79.30010),"4": (19.96839,79.31393),
    "5": (19.96024,79.31746),"6": (19.97239,79.30880),
    "7": (20.00482,79.30233),"8": (19.94894,79.29743),
    "9": (19.93473,79.32349),"10":(19.95375,79.32983),
}
CALIB = {
    "1":{"b0":0.114,"b_raw":1.0026},"2":{"b0":0.231,"b_raw":0.9832},
    "3":{"b0":-0.246,"b_raw":1.0233},"4":{"b0":0.132,"b_raw":0.9851},
    "5":{"b0":0.313,"b_raw":0.9242},"6":{"b0":0.243,"b_raw":0.9593},
    "7":{"b0":0.289,"b_raw":0.9907},"8":{"b0":0.308,"b_raw":1.0085},
    "9":{"b0":0.165,"b_raw":1.0216},"10":{"b0":0.489,"b_raw":0.9847},
}
DEPLOY_START = {
    "1":"2026-08-11","2":"2026-08-11","3":"2026-08-11","4":"2026-08-24",
    "5":"2026-09-04","6":"2026-09-02","7":"2026-08-11","8":"2026-08-11",
    "9":"2026-08-11","10":"2026-08-11",
}
SECTOR_DEG = 45

# ── Stack parameters (CSTPS_dispersion_brief.pdf) ─────────────────────────
Hs = 275.0   # m, physical stack height
Vs = 22.0    # m/s, exit velocity
ds = 6.5     # m, flue diameter
Ts = 400.0   # K, exit gas temperature
Ta = 303.0   # K, ambient

# ── Geometry helpers ───────────────────────────────────────────────────────
def bearing(lat1,lon1,lat2,lon2):
    φ1,φ2=math.radians(lat1),math.radians(lat2)
    Δλ=math.radians(lon2-lon1)
    x=math.sin(Δλ)*math.cos(φ2)
    y=math.cos(φ1)*math.sin(φ2)-math.sin(φ1)*math.cos(φ2)*math.cos(Δλ)
    return (math.degrees(math.atan2(x,y))+360)%360

def haversine_km(lat1,lon1,lat2,lon2):
    R=6371.0; φ1,φ2=math.radians(lat1),math.radians(lat2)
    dφ=math.radians(lat2-lat1); dλ=math.radians(lon2-lon1)
    a=math.sin(dφ/2)**2+math.cos(φ1)*math.cos(φ2)*math.sin(dλ/2)**2
    return R*2*math.asin(math.sqrt(a))

def angle_diff(a,b):
    d=(a-b+360)%360
    return d if d<=180 else d-360

# ── Gaussian dispersion (Pasquill-Gifford Class D, neutral) ───────────────
def sz_pgD(x_km):  return 0.06 * x_km**0.91   # vertical σ, km
def sy_pgD(x_km):  return 0.10 * x_km**0.88   # lateral σ, km

def glc_full(x_arr, He_km):
    """
    Full 2D ground-level centerline concentration (normalized, arbitrary units).
    C(x) ∝ 1/(σy·σz) · exp(-He²/(2·σz²))
    Peaks at shorter distance than crosswind-integrated form because the
    σy factor in the denominator grows faster than σz, pulling the max inward.
    """
    sy = sy_pgD(x_arr)
    sz = sz_pgD(x_arr)
    c  = (1.0 / (sy * sz)) * np.exp(-He_km**2 / (2.0 * sz**2))
    c[x_arr < 0.3] = 0.0
    return c

def plume_rise_briggs(u_ms):
    """Briggs (1975) ultimate buoyant plume rise."""
    rs = ds / 2.0
    Fb = 9.81 * Vs * rs**2 * (Ts - Ta) / Ts
    xf = 120.0 * Fb**0.4 if Fb > 55 else 50.0 * Fb**0.625
    return 1.6 * Fb**(1/3) * xf**(2/3) / u_ms

# ── Load ERA5 ──────────────────────────────────────────────────────────────
print("Loading ERA5…")
era5 = pd.read_csv(ERA5_PATH, parse_dates=["datetime_ist"])
era5 = era5.set_index("datetime_ist").sort_index()
era5 = era5[["wind_dir_850hpa","wind_speed_850hpa_ms"]].resample("h").mean()
mean_u = era5["wind_speed_850hpa_ms"].mean()
print(f"  Mean 850 hPa wind: {mean_u:.1f} m/s")

# ── Load sensors ───────────────────────────────────────────────────────────
def load_sensor(sid):
    df = pd.read_csv(f"{DATA_DIR}/S{sid}_civil-lines.csv")
    ts = next(c for c in df.columns if "time" in c.lower() or "date" in c.lower())
    pm = next(c for c in df.columns if "pm2" in c.lower() or "pm_2" in c.lower())
    df = df.rename(columns={ts:"ts",pm:"pm_raw"})
    df["ts"] = pd.to_datetime(df["ts"],utc=True).dt.tz_localize(None)
    df = df.dropna(subset=["ts","pm_raw"])
    df = df[df["ts"] >= pd.Timestamp(DEPLOY_START[sid])]
    df = df[df["pm_raw"] > 0]
    c = CALIB[sid]; df["pm25"] = (c["b0"] + c["b_raw"] * df["pm_raw"]).clip(lower=0)
    q = df["pm25"].quantile([0.01, 0.99])
    df = df[(df["pm25"] >= q[0.01]) & (df["pm25"] <= q[0.99])]
    return df[["ts","pm25"]].set_index("ts").resample("h").mean()

# ── Wind-sector analysis ───────────────────────────────────────────────────
print("Wind-sector analysis…")
results = []
for sid,(slat,slon) in COORDS.items():
    plat,plon = PLANT
    dist = haversine_km(plat,plon,slat,slon)
    req_wind = (bearing(plat,plon,slat,slon) + 180) % 360

    pm = load_sensor(sid)
    df = pm.join(era5[["wind_dir_850hpa"]], how="inner").dropna()
    df["in_sector"] = df["wind_dir_850hpa"].apply(
        lambda wd: abs(angle_diff(wd,req_wind)) <= SECTOR_DEG/2)
    sect  = df.loc[df["in_sector"],  "pm25"].dropna()
    bkgnd = df.loc[~df["in_sector"], "pm25"].dropna()
    if len(sect)<6 or len(bkgnd)<6: continue

    delta = sect.mean() - bkgnd.mean()
    _,p   = stats.ttest_ind(sect, bkgnd, equal_var=False)
    sem_s = sect.std(ddof=1)/np.sqrt(len(sect))
    sem_b = bkgnd.std(ddof=1)/np.sqrt(len(bkgnd))
    ci95  = 1.96*np.sqrt(sem_s**2+sem_b**2)
    stars = "***" if p<0.001 else ("**" if p<0.01 else ("*" if p<0.05 else "ns"))
    results.append({"station":sid,"dist_km":dist,"delta":delta,
                    "ci95":ci95,"p_val":p,"stars":stars})

df_r = pd.DataFrame(results).sort_values("dist_km").reset_index(drop=True)
print(df_r[["station","dist_km","delta","ci95","stars"]].to_string(index=False))

# ── Fit Gaussian model ─────────────────────────────────────────────────────
# Fit to positive-significant stations (exclude S7 near-field and S5/S10 ns)
pos_sig = df_r[(df_r["delta"] > 0) & (df_r["p_val"] < 0.05)]
x_fit = pos_sig["dist_km"].values
y_fit = pos_sig["delta"].values

def model_func(x, k, He_km):
    return k * glc_full(x, He_km)

popt, _ = optimize.curve_fit(model_func, x_fit, y_fit,
                              p0=[5000, 0.35],
                              bounds=([0, 0.1], [1e7, 0.8]))
k_fit, He_fit = float(popt[0]), float(popt[1])
He_fit_m = He_fit * 1000

print(f"\nFitted He = {He_fit_m:.0f} m  (k = {k_fit:.1f})")
x_model = np.linspace(0.3, 12.0, 500)
C_fit   = model_func(x_model, k_fit, He_fit)
peak_x  = x_model[np.argmax(C_fit)]
print(f"Model peak at {peak_x:.1f} km")

# R² on fit stations
res   = y_fit - model_func(x_fit, k_fit, He_fit)
ss_r  = np.sum(res**2)
ss_t  = np.sum((y_fit - y_fit.mean())**2)
r2    = 1 - ss_r/ss_t
print(f"R² (fit stations) = {r2:.2f}")

# Briggs He for reference
He_briggs = Hs + plume_rise_briggs(mean_u)
He_briggs_km = He_briggs / 1000
C_briggs = model_func(x_model, k_fit * (glc_full(x_fit, He_fit) / glc_full(x_fit, He_briggs_km)).mean(),
                      He_briggs_km)
peak_briggs = x_model[np.argmax(C_briggs)]
print(f"Briggs He = {He_briggs:.0f} m → peak at {peak_briggs:.1f} km")

# Uncertainty band: He ± 60 m around fitted value
He_lo_km = max((He_fit_m - 60)/1000, 0.20)
He_hi_km = (He_fit_m + 60) / 1000
# Rescale band to same peak amplitude as best-fit curve
C_lo = glc_full(x_model.copy(), He_lo_km); C_lo = C_lo / C_lo.max() * C_fit.max()
C_hi = glc_full(x_model.copy(), He_hi_km); C_hi = C_hi / C_hi.max() * C_fit.max()
C_briggs_norm = C_briggs / C_briggs.max() * C_fit.max()

# ── Figure ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9.5, 6.2))
fig.patch.set_facecolor("white")
ax.set_facecolor("#F8F9FA")

ax.axhline(0, color="#444", lw=1.3, zorder=1)

# Near-field zone
ax.axvspan(0, 2.0, alpha=0.13, color="#4FC3F7", zorder=0)
ax.text(1.0, -18.5, "Near-field\n(plume aloft)", ha="center", va="bottom",
        fontsize=8.5, color="#0277BD", style="italic")

# Gaussian model — uncertainty band
ax.fill_between(x_model, C_lo, C_hi, alpha=0.22, color="#E65100", zorder=2,
                label=f"He ± 60 m band ({He_fit_m-60:.0f}–{He_fit_m+60:.0f} m)")
# Best-fit curve
ax.plot(x_model, C_fit, color="#E65100", lw=2.4, ls="-", zorder=3,
        label=f"Gaussian best fit  (He = {He_fit_m:.0f} m, PG-D)")
# Briggs reference
ax.plot(x_model, C_briggs_norm, color="#888", lw=1.5, ls=":", zorder=3,
        label=f"Briggs ΔH (He = {He_briggs:.0f} m)")

# Peak annotation
ax.annotate("", xy=(peak_x, C_fit.max()+1.0), xytext=(peak_x, C_fit.max()+5.5),
            arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.8))
ax.text(peak_x, C_fit.max()+6.2, f"Peak ≈ {peak_x:.1f} km\n(brief: 4.8 km)",
        ha="center", va="bottom", fontsize=8.5, color="#E65100", fontweight="bold")

# Observed data points
c_pos_sig = "#C62828"
c_low_sig = "#1565C0"   # S7
c_ns      = "#757575"

for _, row in df_r.iterrows():
    sid = row["station"]
    x, y, ci = row["dist_km"], row["delta"], row["ci95"]
    sig, p = row["stars"], row["p_val"]
    is_neg_nf = (y < 0 and x < 2.0)
    clr = c_low_sig if is_neg_nf else (c_pos_sig if sig != "ns" else c_ns)
    filled = (sig != "ns") or is_neg_nf

    ax.errorbar(x, y, yerr=ci, fmt="none", ecolor=clr,
                elinewidth=1.5, capsize=4.5, capthick=1.5, zorder=5)
    ax.scatter(x, y, s=95,
               color=clr if filled else "white",
               edgecolors=clr, linewidths=2.0, zorder=6)

    # station label
    ha = "right" if x < 5.5 else "left"
    dx = -0.18 if ha == "right" else 0.18
    dy = 1.2 if y > 0 else -2.2
    ax.annotate(f"S{sid}", (x+dx, y+dy), fontsize=8.5, fontweight="bold",
                color=clr, ha=ha,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.7, ec="none"))

    # significance stars
    if sig != "ns":
        ax.text(x, y+ci+0.8, sig, ha="center", va="bottom",
                fontsize=7.5, color=clr, fontweight="bold")

# Axes
ax.set_xlabel("Distance from CSTPS (km)", fontsize=12)
ax.set_ylabel(r"$\Delta$PM$_{2.5}$ (µg/m³)   [plant-sector − background]", fontsize=11)
ax.set_title(
    "Ground-Level PM$_{2.5}$ Enhancement vs. Distance from CSTPS\n"
    "Observed Wind-Sector Differences vs. Gaussian Plume Prediction  (PG Class D, full centerline)",
    fontsize=11.5, fontweight="bold", pad=12)

ax.set_xlim(0.0, 11.0)
ypad = 4.0
ax.set_ylim(-22, max(df_r["delta"].max()+df_r["ci95"].max()+ypad, C_fit.max()+10))
ax.set_xticks(range(0, 12))
ax.grid(axis="y", ls="--", lw=0.6, alpha=0.45, color="#ccc")
ax.grid(axis="x", ls=":", lw=0.4, alpha=0.3, color="#ccc")

# Legend
legend_handles = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=c_pos_sig,
           markersize=9, markeredgewidth=1.8, label='Significant positive Δ (p < 0.05)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='white',
           markersize=9, markeredgewidth=2.0, markeredgecolor=c_ns,
           label='Non-significant (p ≥ 0.05)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=c_low_sig,
           markersize=9, markeredgewidth=1.8, label='Near-field suppression (S7, 1.1 km)'),
    Line2D([0],[0], ls='-', color='#E65100', lw=2.4,
           label=f'Gaussian best fit (He = {He_fit_m:.0f} m, PG-D)'),
    mpatches.Patch(facecolor='#E65100', alpha=0.25,
                   label=f'He ± 60 m uncertainty band'),
    Line2D([0],[0], ls=':', color='#888', lw=1.6,
           label=f'Briggs ΔH estimate (He = {He_briggs:.0f} m)'),
]
ax.legend(handles=legend_handles, loc="upper right", fontsize=8.5,
          framealpha=0.92, edgecolor="#ccc", handlelength=2.0)

# Caption
dH_briggs = He_briggs - Hs
caption = (
    f"Stack: H_s = {Hs:.0f} m, V_s = {Vs:.0f} m/s, d = {ds:.1f} m, T_s = {Ts:.0f} K  |  "
    f"Briggs ΔH = {dH_briggs:.0f} m at ū = {mean_u:.1f} m/s  |  "
    f"Fitted He = {He_fit_m:.0f} m (peak at {peak_x:.1f} km; brief ensemble: 4.8 km)  |  "
    f"10 LCS stations, Aug–Oct 2026"
)
fig.text(0.5, -0.03, caption, ha="center", fontsize=7.5, color="#555")

plt.tight_layout(rect=[0, 0, 1, 1])
plt.savefig(FIG_OUT, dpi=200, bbox_inches="tight")
print(f"\nSaved → {FIG_OUT}")
