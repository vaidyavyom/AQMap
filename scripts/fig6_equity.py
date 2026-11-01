#!/usr/bin/env python3
"""
fig6_equity.py  —  Environmental equity figure (single panel).
Shows that the highest-PM₂.₅ stations sit in Chandrapur's lower-income
peri-urban fringe, while the established civil-lines/urban-core area is cleaner.

Output: figures/fig6_equity.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import os

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data", "sensors")
FIG_OUT  = os.path.join(ROOT, "figures", "fig6_equity.png")

# ── Sensor metadata ────────────────────────────────────────────────────────
PLANT = (20.00837, 79.29206)

COORDS = {
    "1":(19.98062,79.23735),"2":(19.97237,79.28563),"3":(19.98250,79.30010),
    "4":(19.96839,79.31393),"5":(19.96024,79.31746),"6":(19.97239,79.30880),
    "7":(20.00482,79.30233),"8":(19.94894,79.29743),"9":(19.93473,79.32349),
    "10":(19.95375,79.32983),
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
CSV_MAP = {k: f"S{k}_civil-lines.csv" for k in COORDS}
MONSOON_END = pd.Timestamp("2026-09-22")

# Neighbourhood characterisation — based on Chandrapur urban geography
# (civil lines = established admin / professional quarter;
#  peri-urban SE = MIDC-adjacent labour colonies & migrant settlements)
NEIGHBOURHOOD = {
    "1":  "Civil Lines\n(established, professional)",
    "2":  "Urban core\n(commercial / residential)",
    "3":  "Urban core\n(commercial / residential)",
    "4":  "SE peri-urban\n(labour colony / informal)",
    "5":  "SE peri-urban\n(labour colony / informal)",
    "6":  "Urban core\n(commercial / residential)",
    "7":  "Industrial buffer\n(near CSTPS fence)",
    "8":  "SE peri-urban\n(labour colony / informal)",
    "9":  "SE fringe\n(peri-urban / informal)",
    "10": "SE peri-urban\n(labour colony / informal)",
}

# Group colours — three tiers of income
TIER_COLOR = {
    "Civil Lines\n(established, professional)":  "#1565C0",   # dark blue
    "Urban core\n(commercial / residential)":    "#2E7D32",   # dark green
    "Industrial buffer\n(near CSTPS fence)":     "#E65100",   # deep orange
    "SE peri-urban\n(labour colony / informal)": "#B71C1C",   # deep red
    "SE fringe\n(peri-urban / informal)":        "#C62828",   # red
}

# ── Load sensor data ───────────────────────────────────────────────────────
def load(sid):
    df = pd.read_csv(f"{DATA_DIR}/{CSV_MAP[sid]}")
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
STA = {}
for sid in COORDS:
    df   = load(sid)
    post = df[df.index > MONSOON_END]["pm25"]
    STA[sid] = {
        "pm":  post.mean(),
        "sem": post.sem(),
        "n":   len(post),
        "lat": COORDS[sid][0],
        "lon": COORDS[sid][1],
        "nb":  NEIGHBOURHOOD[sid],
        "col": TIER_COLOR[NEIGHBOURHOOD[sid]],
    }


# ── Figure: 1×2 layout  (map | ranked comparison) ─────────────────────────
vmin, vmax = 38, 56

print("Building figure…")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.linewidth": 0.8,
})

fig, ax_bar = plt.subplots(1, 1, figsize=(9, 6.5))
fig.patch.set_facecolor("white")

# ══════════════════════════════════════════════════════════════════════════
# RIGHT: Ranked horizontal bar — sorted by PM2.5, coloured by tier
# ══════════════════════════════════════════════════════════════════════════
ax = ax_bar

# Sort stations highest→lowest PM2.5
order = sorted(STA.keys(), key=lambda s: STA[s]["pm"], reverse=True)

y   = np.arange(len(order))
pms = np.array([STA[s]["pm"]  for s in order])
sem = np.array([STA[s]["sem"] for s in order])
col = [STA[s]["col"] for s in order]

# Horizontal bars
ci95 = sem * 1.96
bars = ax.barh(y, pms, xerr=ci95, color=col, alpha=0.88,
               edgecolor="white", linewidth=0.8, height=0.65,
               error_kw={"ecolor":"#222","lw":1.3,"capsize":4}, zorder=3)

# WHO / NAAQS reference lines (clipped to axes)
ax.axvline(60, color="#B71C1C", ls="--", lw=1.2, alpha=0.7, zorder=2, clip_on=True)

# Station labels on left
ax.set_yticks(y)
labels = [f"S{s}" for s in order]
ax.set_yticklabels(labels, fontsize=10, fontweight="bold")

# Neighbourhood annotation on the right of each bar
for i, sid in enumerate(order):
    nb   = STA[sid]["nb"]
    pm_v = STA[sid]["pm"]
    # short neighbourhood name (first line only)
    short = nb.split("\n")[0]
    ax.text(pm_v + ci95[i] + 0.5, i, short,
            va="center", ha="left", fontsize=8.2, color=col[i])

# Reference line at Civil Lines / core average
avg_urban = np.mean([STA[s]["pm"] for s in ["1","2","3"]])
ax.axvline(avg_urban, color="#1565C0", ls="-.", lw=1.5, alpha=0.7, zorder=2)
ax.text(avg_urban - 0.3, 0.3,
        f"Established urban avg\n{avg_urban:.1f} µg/m³",
        ha="right", va="top", fontsize=8, color="#1565C0", style="italic")

# NAAQS label inside axes
ax.text(60.3, len(order) - 0.8, "NAAQS\n60 µg", fontsize=7.5,
        color="#B71C1C", va="top", ha="left", rotation=0)

ax.set_xlabel("Post-monsoon mean PM₂.₅ (µg/m³)  ±  95% CI", fontsize=10)
ax.set_xlim(vmin - 7, 67)
ax.set_ylim(-0.5, len(order) - 0.5)
ax.invert_yaxis()
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.grid(axis="x", color="#ddd", lw=0.7, zorder=0)

# Equity gap stat box
avg_peri  = np.mean([STA[s]["pm"] for s in ["4","5","8","9","10"]])
avg_urban = np.mean([STA[s]["pm"] for s in ["1","2","3"]])
ax.text(0.98, 0.03,
        f"SE peri-urban avg:  {avg_peri:.1f} µg/m³\n"
        f"Civil Lines/core avg:  {avg_urban:.1f} µg/m³\n"
        f"Difference:  +{avg_peri-avg_urban:.1f} µg/m³  ({(avg_peri/avg_urban-1)*100:.0f}% higher)",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=9, color="#333",
        bbox=dict(boxstyle="round,pad=0.5", fc="#FFF9C4", ec="#F9A825", lw=1.2))

ax.set_title("Stations ranked by post-monsoon PM₂.₅ exposure",
             fontsize=11, fontweight="bold", loc="left", pad=8)

# ── Shared suptitle ────────────────────────────────────────────────────────
fig.suptitle(
    "Environmental Equity  —  Chandrapur's lower-income peri-urban fringe "
    "bears the highest PM₂.₅ burden",
    fontsize=12.5, fontweight="bold", y=1.00,
)

plt.savefig(FIG_OUT, dpi=220, bbox_inches="tight")
print(f"Saved → {FIG_OUT}")
