#!/usr/bin/env python3
"""
fig1_spatial_maps.py  —  Proportional-symbol spatial map (no IDW).
Shows actual station PM2.5 values as sized, colour-coded bubbles on a
basemap.  No spatial interpolation — with 10 stations IDW produces
meaningless colour gradients.

Neighbourhood colours match fig6_equity.py so the two figures
cross-reference visually.

Run from AQMap root:  python3 scripts/fig1_spatial_maps.py
Output:               figures/fig1_spatial_maps.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import os

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data", "sensors")
FIG_OUT  = os.path.join(ROOT, "figures", "fig1_spatial_maps.png")

# ── Coordinates ────────────────────────────────────────────────────────────
PLANT = (20.00837, 79.29206)   # CSTPS centroid

COORDS = {
    "1": (19.98062, 79.23735), "2": (19.97237, 79.28563),
    "3": (19.98250, 79.30010), "4": (19.96839, 79.31393),
    "5": (19.96024, 79.31746), "6": (19.97239, 79.30880),
    "7": (20.00482, 79.30233), "8": (19.94894, 79.29743),
    "9": (19.93473, 79.32349), "10":(19.95375, 79.32983),
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

MONSOON_END = pd.Timestamp("2026-09-22")
lat_min, lat_max = 19.91, 20.03
lon_min, lon_max = 79.20, 79.37

# ── Neighbourhood colours — identical to fig6_equity.py ───────────────────
NEIGHBOURHOOD = {
    "1":  "Civil Lines",
    "2":  "Urban core",
    "3":  "Urban core",
    "4":  "SE peri-urban",
    "5":  "SE peri-urban",
    "6":  "Urban core",
    "7":  "Industrial buffer",
    "8":  "SE peri-urban",
    "9":  "SE fringe",
    "10": "SE peri-urban",
}

NB_COLOR = {
    "Civil Lines":        "#1565C0",   # dark blue
    "Urban core":         "#2E7D32",   # dark green
    "Industrial buffer":  "#E65100",   # deep orange
    "SE peri-urban":      "#B71C1C",   # deep red
    "SE fringe":          "#C62828",   # red
}

# ── Load sensor data ────────────────────────────────────────────────────────
def load(sid):
    df = pd.read_csv(f"{DATA_DIR}/S{sid}_civil-lines.csv")
    ts = next(c for c in df.columns if "time" in c.lower() or "date" in c.lower())
    pm = next(c for c in df.columns if "pm2" in c.lower() or "pm_2" in c.lower())
    df = df.rename(columns={ts:"ts", pm:"pm_raw"})
    df["ts"] = pd.to_datetime(df["ts"], utc=True).dt.tz_localize(None)
    df = df.dropna(subset=["ts","pm_raw"])
    df = df[df["ts"] >= pd.Timestamp(DEPLOY_START[sid])]
    df = df[df["pm_raw"] > 0]
    c  = CALIB[sid]
    df["pm25"] = (c["b0"] + c["b_raw"] * df["pm_raw"]).clip(lower=0)
    q = df["pm25"].quantile([0.01, 0.99])
    return df[(df["pm25"] >= q[0.01]) & (df["pm25"] <= q[0.99])].set_index("ts")

print("Loading sensors…")
sm, sp = {}, {}
for sid in COORDS:
    df = load(sid)
    sm[sid] = df[df.index <= MONSOON_END]["pm25"].mean()
    post     = df[df.index >  MONSOON_END]["pm25"]
    sp[sid]  = post.mean() if len(post) > 0 else np.nan
    print(f"  S{sid:>2}: monsoon={sm[sid]:.1f}  post-monsoon={sp[sid]:.1f} µg/m³")

# ── Try satellite basemap (needs contextily + internet) ─────────────────────
USE_SATELLITE = False
try:
    import contextily as ctx
    USE_SATELLITE = True
    print("contextily available — will fetch satellite tiles.")
except ImportError:
    print("contextily not installed — using plain basemap.")

# ── Figure setup ────────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.linewidth": 0.8,
})

fig, axes = plt.subplots(1, 2, figsize=(14, 6.8))
fig.patch.set_facecolor("#f8f8f8" if not USE_SATELLITE else "white")

outline_w = [pe.withStroke(linewidth=2.5, foreground="black")]
outline_b = [pe.withStroke(linewidth=2.5, foreground="white")]

# Bubble size scale: area ∝ PM2.5 value
# Reference: a circle of 900 pt² represents 45 µg/m³
REF_VAL  = 45.0
REF_SIZE = 900

def pm_to_size(val):
    """Circle area proportional to PM2.5."""
    return REF_SIZE * (val / REF_VAL) ** 1.5

panels = [
    (axes[0], sm,
     "(a)  Monsoon Period\n(Aug 11 – Sep 22, 2026)"),
    (axes[1], sp,
     "(b)  Post-Monsoon Period\n(Sep 23 – Oct 6, 2026)"),
]

for ax, svals, title in panels:
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    # ── basemap ─────────────────────────────────────────────────────────────
    if USE_SATELLITE:
        try:
            ctx.add_basemap(ax, crs="EPSG:4326",
                            source=ctx.providers.Esri.WorldImagery,
                            zoom=13, attribution=False)
            bg_dark = True
        except Exception as e:
            print(f"  Tile fetch failed ({e}), using plain basemap.")
            ax.set_facecolor("#e8ecf0")
            bg_dark = False
    else:
        ax.set_facecolor("#e8ecf0")   # light blue-grey geographic feel
        bg_dark = False

    # ── light reference grid ─────────────────────────────────────────────────
    if not bg_dark:
        for glon in np.arange(79.20, 79.38, 0.03):
            ax.axvline(glon, color="white", lw=0.5, alpha=0.5, zorder=1)
        for glat in np.arange(19.91, 20.04, 0.02):
            ax.axhline(glat, color="white", lw=0.5, alpha=0.5, zorder=1)

    outline = outline_b if bg_dark else outline_w

    # ── station bubbles ──────────────────────────────────────────────────────
    # Draw in two passes: shadow ring first, then coloured fill
    for sid, (slat, slon) in COORDS.items():
        val = svals.get(sid, np.nan)
        if np.isnan(val):
            continue
        col  = NB_COLOR[NEIGHBOURHOOD[sid]]
        sz   = pm_to_size(val)

        # Shadow ring
        ax.scatter(slon, slat, s=sz * 1.25,
                   color="black", alpha=0.25, zorder=4, linewidths=0)
        # Filled bubble
        ax.scatter(slon, slat, s=sz,
                   color=col, edgecolors="white", linewidths=1.5,
                   alpha=0.90, zorder=5)
        # Station label + value
        label = f"S{sid}\n{val:.1f} µg/m³"
        va_off = 0.006 if slat < PLANT[0] else -0.006
        va_align = "bottom" if va_off > 0 else "top"
        ax.annotate(
            label,
            (slon, slat + va_off),
            ha="center", va=va_align,
            fontsize=7.5, fontweight="bold",
            color="white" if bg_dark else "#111",
            path_effects=outline if bg_dark else [],
            zorder=6,
        )

    # ── CSTPS plant marker ───────────────────────────────────────────────────
    ax.scatter(PLANT[1], PLANT[0],
               marker="s", s=180,
               color="#C62828", edgecolors="white", linewidths=2.0, zorder=8)
    ax.annotate(
        "CSTPS\n(power plant)",
        (PLANT[1], PLANT[0]),
        xytext=(0.006, 0.006), textcoords="offset points",
        ha="left", va="bottom",
        fontsize=8.5, fontweight="bold",
        color="white" if bg_dark else "#C62828",
        path_effects=outline_b if bg_dark else [],
        zorder=8,
    )

    # ── Plume transport arrow (WNW→ESE dominant direction) ───────────────────
    ax.annotate(
        "",
        xy=(PLANT[1] + 0.030, PLANT[0] - 0.022),
        xytext=(PLANT[1] - 0.010, PLANT[0] + 0.008),
        arrowprops=dict(
            arrowstyle="-|>",
            color="white" if bg_dark else "#555",
            lw=2.0, mutation_scale=16,
        ),
        zorder=7,
    )
    ax.text(
        PLANT[1] + 0.033, PLANT[0] - 0.025,
        "Prevailing\nplume axis",
        fontsize=7.5,
        color="white" if bg_dark else "#444",
        ha="left", va="top",
        path_effects=outline_b if bg_dark else [],
        zorder=7,
    )

    # ── Network mean box ─────────────────────────────────────────────────────
    net_mean = np.nanmean(list(svals.values()))
    ax.text(
        0.02, 0.03, f"Network mean: {net_mean:.1f} µg/m³",
        transform=ax.transAxes,
        fontsize=9, fontweight="bold",
        color="white" if bg_dark else "#222",
        bbox=dict(boxstyle="round,pad=0.4",
                  fc="black" if bg_dark else "white",
                  alpha=0.65, ec="none"),
        zorder=9,
    )

    # ── NAAQS reference note ─────────────────────────────────────────────────
    ax.text(
        0.98, 0.03, "NAAQS 24h = 60 µg/m³",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7.5, style="italic",
        color="#bbb" if bg_dark else "#888",
        zorder=9,
    )

    ax.set_title(title, fontsize=11, fontweight="bold", pad=8, loc="left")
    ax.set_xlabel("Longitude (°E)", fontsize=9)
    ax.set_ylabel("Latitude (°N)", fontsize=9)
    ax.tick_params(labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor("#aaa")

# ── Legend ───────────────────────────────────────────────────────────────────
legend_entries = [
    mpatches.Patch(color="#C62828", label="CSTPS — coal power plant"),
    mpatches.Patch(color=NB_COLOR["Civil Lines"],       label="Civil Lines (established, professional)"),
    mpatches.Patch(color=NB_COLOR["Urban core"],        label="Urban core (commercial / residential)"),
    mpatches.Patch(color=NB_COLOR["Industrial buffer"], label="Industrial buffer (near CSTPS fence)"),
    mpatches.Patch(color=NB_COLOR["SE peri-urban"],     label="SE peri-urban (labour colony / informal)"),
    mpatches.Patch(color=NB_COLOR["SE fringe"],         label="SE fringe (peri-urban / informal)"),
]

# Bubble size reference
for ref_val in [35, 45, 55]:
    sz = pm_to_size(ref_val)
    legend_entries.append(
        Line2D([0],[0], marker="o", color="none",
               markerfacecolor="#888", markeredgecolor="#555",
               markersize=np.sqrt(sz) * 0.45,
               label=f"{ref_val} µg/m³")
    )

fig.legend(handles=legend_entries, loc="lower center", ncol=5,
           fontsize=8, framealpha=0.92, edgecolor="#ccc",
           bbox_to_anchor=(0.5, -0.01))

# ── Note on method ─────────────────────────────────────────────────────────
fig.text(
    0.5, -0.055,
    "Circle area ∝ mean PM₂.₅  ·  No spatial interpolation applied"
    "  ·  Colour indicates neighbourhood classification",
    ha="center", va="top", fontsize=8, color="#555", style="italic",
)

fig.suptitle(
    "PM₂.₅ Spatial Distribution — Chandrapur, Maharashtra, 2026\n"
    "10-Station Low-Cost Sensor Network",
    fontsize=13, fontweight="bold", y=1.01,
)

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig(FIG_OUT, dpi=220, bbox_inches="tight")
print(f"\nSaved → {FIG_OUT}")
