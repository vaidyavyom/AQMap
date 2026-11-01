#!/usr/bin/env python3
"""
fig_rain_rebound.py  —  Rain rebound / wet-scavenging event analysis.

Identifies rain events from the sensor network itself: wet scavenging
produces rapid, coherent PM2.5 drops across ≥6 stations simultaneously.
For each event, tracks PM2.5 in 24 h before, during, and 48 h after.
Tests whether post-rain PM2.5 rebounds above the pre-rain baseline
(the "pollution rebound paradox" documented in Delhi monsoon 2026).

Run from AQMap root:  python3 scripts/fig_rain_rebound.py
Outputs:
    figures/fig_rain_rebound.png
    data/rain_events.csv
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import os

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data", "sensors")
FIG_OUT  = os.path.join(ROOT, "figures", "fig_rain_rebound.png")
CSV_OUT  = os.path.join(ROOT, "data", "rain_events.csv")

# ── Sensor metadata ────────────────────────────────────────────────────────
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
MONSOON_END = pd.Timestamp("2026-09-22")

# ── Load sensors ─────────────────────────────────────────────────────────
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
    df = df[(df["pm25"] >= q[0.01]) & (df["pm25"] <= q[0.99])]
    return df[["ts","pm25"]].set_index("ts").resample("1h").mean()

print("Loading sensors…")
station_dfs = {}
for sid in COORDS:
    station_dfs[sid] = load(sid)
    print(f"  S{sid}: {len(station_dfs[sid])} hourly records")

# Network mean (mean across available stations each hour)
combined = pd.concat([df.rename(columns={"pm25": sid}) for sid, df in station_dfs.items()], axis=1)
network  = combined.mean(axis=1).rename("pm25_net")

# Station availability per hour
n_avail  = combined.notna().sum(axis=1)

# Restrict to monsoon period for event detection (rain events are monsoon)
mon_net = network[network.index <= MONSOON_END].copy()

print(f"\nNetwork mean during monsoon: {mon_net.mean():.1f} µg/m³")

# ── Rain event detection ─────────────────────────────────────────────────
# A wet-scavenging event shows:
#   1. Rapid drop: >20% fall in rolling 3h mean, absolute drop >5 µg/m³
#   2. Coherent: drop simultaneous across ≥5 stations
#   3. Duration: lasts ≥2 hours

ROLL_WIN   = 3    # hours for rolling mean
DROP_PCT   = 0.20 # fractional drop threshold
DROP_ABS   = 5.0  # µg/m³ minimum absolute drop
MIN_STNS   = 5    # minimum stations showing simultaneous drop
MIN_HOURS  = 2    # minimum event duration

roll = mon_net.rolling(ROLL_WIN, center=True, min_periods=2).mean()
# Change over previous 3 hours
delta_3h = roll.diff(3)           # absolute change
frac_3h  = delta_3h / roll.shift(3)  # fractional change

# Flag hours that look like active wet scavenging
scav_flag = (frac_3h < -DROP_PCT) & (delta_3h < -DROP_ABS)

# Confirm with station-level coherence
def count_dropping_stations(t):
    row = combined.loc[t] if t in combined.index else pd.Series(dtype=float)
    if len(row) < MIN_STNS:
        return 0
    prev = combined.loc[:t].tail(3).mean()
    drop = (row - prev) / prev.where(prev > 0)
    return (drop < -0.15).sum()

print("Checking station coherence for flagged hours…")
confirmed = []
for t in scav_flag[scav_flag].index:
    n_drop = count_dropping_stations(t)
    if n_drop >= MIN_STNS:
        confirmed.append(t)

print(f"  Confirmed scavenging hours: {len(confirmed)}")

# Group confirmed hours into events (gap < 3h → same event)
events = []
if confirmed:
    cur_start = confirmed[0]
    cur_end   = confirmed[0]
    for t in confirmed[1:]:
        if (t - cur_end).total_seconds() / 3600 <= 3:
            cur_end = t
        else:
            if (cur_end - cur_start).total_seconds() / 3600 >= MIN_HOURS - 1:
                events.append((cur_start, cur_end))
            cur_start = t
            cur_end   = t
    if (cur_end - cur_start).total_seconds() / 3600 >= MIN_HOURS - 1:
        events.append((cur_start, cur_end))

print(f"  Rain events identified: {len(events)}")
for i, (s, e) in enumerate(events):
    dur = (e - s).total_seconds() / 3600 + 1
    pm_during = network.loc[s:e].mean()
    print(f"    Event {i+1}: {s.strftime('%b %d %H:%M')} → {e.strftime('%b %d %H:%M')} ({dur:.0f}h), PM={pm_during:.1f} µg/m³")

# ── Per-event composite windows ───────────────────────────────────────────
PRE_H  = 24   # hours before event start
POST_H = 48   # hours after event end

event_records = []
composites    = []   # aligned time series, t=0 = event start

for i, (ev_start, ev_end) in enumerate(events):
    dur_h = int((ev_end - ev_start).total_seconds() / 3600) + 1

    pre_start  = ev_start - pd.Timedelta(hours=PRE_H)
    post_end   = ev_end   + pd.Timedelta(hours=POST_H)

    window = network.loc[pre_start:post_end]
    if len(window) < PRE_H + dur_h + POST_H // 2:
        print(f"  Event {i+1}: insufficient data in window, skipping")
        continue

    pre_mean    = network.loc[pre_start : ev_start - pd.Timedelta(hours=1)].mean()
    during_mean = network.loc[ev_start  : ev_end].mean()
    post_24h    = network.loc[ev_end + pd.Timedelta(hours=1) : ev_end + pd.Timedelta(hours=24)].mean()
    post_48h    = network.loc[ev_end + pd.Timedelta(hours=1) : ev_end + pd.Timedelta(hours=48)].mean()

    drop_abs  = during_mean - pre_mean
    drop_pct  = drop_abs / pre_mean * 100 if pre_mean > 0 else np.nan
    reb_24_abs = post_24h - pre_mean
    reb_24_pct = reb_24_abs / pre_mean * 100 if pre_mean > 0 else np.nan
    reb_48_abs = post_48h - pre_mean
    reb_48_pct = reb_48_abs / pre_mean * 100 if pre_mean > 0 else np.nan

    event_records.append({
        "event": i+1,
        "start": ev_start, "end": ev_end, "duration_h": dur_h,
        "pre_mean": pre_mean, "during_mean": during_mean,
        "post_24h_mean": post_24h, "post_48h_mean": post_48h,
        "drop_abs": drop_abs, "drop_pct": drop_pct,
        "rebound_24h_abs": reb_24_abs, "rebound_24h_pct": reb_24_pct,
        "rebound_48h_abs": reb_48_abs, "rebound_48h_pct": reb_48_pct,
    })

    # Aligned time series (hourly offset from event start)
    aligned = []
    for h in range(-PRE_H, POST_H + dur_h + 1):
        t = ev_start + pd.Timedelta(hours=h)
        v = network.get(t, np.nan)
        aligned.append({"hour": h, "pm25": v, "norm": v / pre_mean if pre_mean > 0 else np.nan})
    composites.append(pd.DataFrame(aligned).set_index("hour"))

df_events = pd.DataFrame(event_records)
df_events.to_csv(CSV_OUT, index=False)
print(f"\nEvent table saved → {CSV_OUT}")
print(df_events[["event","start","pre_mean","during_mean","post_24h_mean","post_48h_mean","drop_pct","rebound_24h_pct","rebound_48h_pct"]].to_string())

# ── Composite mean across events ──────────────────────────────────────────
if composites:
    comp_all = pd.concat(composites, axis=0).reset_index()
    comp_mean = comp_all.groupby("hour")["norm"].agg(["mean","sem"]).reset_index()
    comp_pm   = comp_all.groupby("hour")["pm25"].agg(["mean","sem"]).reset_index()
else:
    comp_mean = pd.DataFrame()
    print("No valid events for composite — check data coverage.")

# ── Figure ────────────────────────────────────────────────────────────────
matplotlib.rcParams.update({"font.family":"DejaVu Sans","font.size":9.5,"axes.linewidth":0.8})

n_events = len(df_events)
has_composite = len(comp_mean) > 0

fig = plt.figure(figsize=(14, 9))
fig.patch.set_facecolor("white")

gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.30)
ax_ts   = fig.add_subplot(gs[0, :])   # top full-width: full monsoon timeseries
ax_comp = fig.add_subplot(gs[1, 0])   # bottom-left: composite event profile
ax_bar  = fig.add_subplot(gs[1, 1])   # bottom-right: per-event bar summary

# ── Panel A: full monsoon timeseries with event shading ────────────────────
ax_ts.plot(mon_net.index, mon_net.values, color="#2E7D32", lw=1.0, alpha=0.7, label="Network mean PM₂.₅")
roll_plot = mon_net.rolling(6, center=True, min_periods=3).mean()
ax_ts.plot(roll_plot.index, roll_plot.values, color="#1B5E20", lw=1.8, label="6-h rolling mean")

event_colors = plt.cm.Set2.colors
for i, (ev_start, ev_end) in enumerate(events):
    col = event_colors[i % len(event_colors)]
    ax_ts.axvspan(ev_start, ev_end, alpha=0.25, color=col, zorder=0)
    ax_ts.axvspan(ev_end, ev_end + pd.Timedelta(hours=48), alpha=0.10, color=col, zorder=0)
    ax_ts.annotate(f"E{i+1}", xy=(ev_start, mon_net.max() * 0.95),
                   fontsize=8, color=col, fontweight="bold", ha="left")

ax_ts.axhline(15, color="#B71C1C", ls="--", lw=1.0, alpha=0.6, label="WHO 15 µg/m³")
ax_ts.set_ylabel("Network mean PM₂.₅ (µg/m³)", fontsize=9)
ax_ts.set_title("(a)  Monsoon-period PM₂.₅ with identified wet-scavenging events (shaded) and 48-h post-rain window",
                fontsize=10, fontweight="bold", loc="left", pad=6)
ax_ts.legend(fontsize=8, loc="upper right", framealpha=0.85)
ax_ts.spines["top"].set_visible(False)
ax_ts.spines["right"].set_visible(False)
ax_ts.tick_params(axis="x", labelsize=8)
ax_ts.tick_params(axis="y", labelsize=8)

# ── Panel B: composite aligned event profile ────────────────────────────
if has_composite:
    ev_dur_med = int(df_events["duration_h"].median())

    ax_comp.axvspan(-PRE_H, 0, alpha=0.06, color="#1565C0", label="Pre-rain window")
    ax_comp.axvspan(0, ev_dur_med, alpha=0.12, color="#42A5F5", label="Median event")
    ax_comp.axvspan(ev_dur_med, ev_dur_med + POST_H, alpha=0.06, color="#E65100", label="Post-rain window")

    ax_comp.axhline(1.0, color="#1565C0", ls="-.", lw=1.2, alpha=0.7, label="Pre-rain baseline")

    # Individual event traces
    for comp in composites:
        ax_comp.plot(comp.index, comp["norm"].values, color="#aaa", lw=0.8, alpha=0.5)

    # Composite mean ± SE
    ax_comp.fill_between(comp_mean["hour"],
                         comp_mean["mean"] - comp_mean["sem"],
                         comp_mean["mean"] + comp_mean["sem"],
                         alpha=0.25, color="#E65100")
    ax_comp.plot(comp_mean["hour"], comp_mean["mean"],
                 color="#E65100", lw=2.2, label="Composite mean ± SE")

    # Mark key checkpoints
    ax_comp.axvline(0, color="#42A5F5", lw=1.4, ls="--", alpha=0.8)
    ax_comp.axvline(ev_dur_med, color="#1565C0", lw=1.4, ls="--", alpha=0.8)
    ax_comp.axvline(ev_dur_med + 24, color="#E65100", lw=1.0, ls=":", alpha=0.7)
    ax_comp.axvline(ev_dur_med + 48, color="#B71C1C", lw=1.0, ls=":", alpha=0.7)

    ax_comp.set_xlabel("Hours relative to event start", fontsize=9)
    ax_comp.set_ylabel("PM₂.₅ normalised to pre-rain mean", fontsize=9)
    ax_comp.set_title("(b)  Composite event profile (normalised)", fontsize=10, fontweight="bold", loc="left", pad=6)
    ax_comp.legend(fontsize=7.5, loc="upper right", framealpha=0.85)
    ax_comp.set_xlim(-PRE_H, ev_dur_med + POST_H)
    ax_comp.spines["top"].set_visible(False)
    ax_comp.spines["right"].set_visible(False)
else:
    ax_comp.text(0.5, 0.5, "Insufficient events\nfor composite", ha="center", va="center",
                 transform=ax_comp.transAxes, fontsize=11, color="#888")

# ── Panel C: per-event summary bars ────────────────────────────────────────
if n_events > 0:
    x = np.arange(n_events)
    w = 0.25
    pre  = df_events["pre_mean"].values
    dur  = df_events["during_mean"].values
    p24  = df_events["post_24h_mean"].values
    p48  = df_events["post_48h_mean"].values

    ax_bar.bar(x - w,   pre, w, label="Pre-rain (24h mean)",    color="#1565C0", alpha=0.85)
    ax_bar.bar(x,       dur, w, label="During rain",             color="#42A5F5", alpha=0.85)
    ax_bar.bar(x + w,   p24, w, label="Post-rain 24h mean",      color="#FF7043", alpha=0.85)
    ax_bar.bar(x + 2*w, p48, w, label="Post-rain 48h mean",      color="#B71C1C", alpha=0.85)

    # Rebound annotation: arrow if post-48h > pre
    for i, row in df_events.iterrows():
        if row["rebound_48h_pct"] > 0:
            ax_bar.annotate(f"+{row['rebound_48h_pct']:.0f}%",
                           xy=(i + 2*w, row["post_48h_mean"] + 0.5),
                           ha="center", va="bottom", fontsize=7.5,
                           color="#B71C1C", fontweight="bold")
        else:
            ax_bar.annotate(f"{row['rebound_48h_pct']:.0f}%",
                           xy=(i + 2*w, row["post_48h_mean"] + 0.5),
                           ha="center", va="bottom", fontsize=7.5,
                           color="#2E7D32")

    ax_bar.set_xticks(x + w/2)
    ax_bar.set_xticklabels([f"Event {int(r['event'])}\n({r['start'].strftime('%b %d')})"
                             for _, r in df_events.iterrows()], fontsize=8)
    ax_bar.set_ylabel("PM₂.₅ (µg/m³)", fontsize=9)
    ax_bar.set_title("(c)  Per-event PM₂.₅: pre / during / post-rain", fontsize=10, fontweight="bold", loc="left", pad=6)
    ax_bar.legend(fontsize=7.5, loc="upper right", framealpha=0.85)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.axhline(15, color="#B71C1C", ls="--", lw=0.8, alpha=0.5)

    # Summary stat box
    mean_reb_48 = df_events["rebound_48h_pct"].mean()
    n_above_pre = (df_events["rebound_48h_pct"] > 0).sum()
    ax_bar.text(0.02, 0.97,
                f"Events with post-48h rebound\nabove pre-rain baseline: {n_above_pre}/{n_events}\n"
                f"Mean 48h rebound: {mean_reb_48:+.1f}%",
                transform=ax_bar.transAxes, va="top", ha="left",
                fontsize=8.5, color="#333",
                bbox=dict(boxstyle="round,pad=0.4", fc="#FFF9C4", ec="#F9A825", lw=1.1))
else:
    ax_bar.text(0.5, 0.5, "No events detected", ha="center", va="center",
                transform=ax_bar.transAxes, fontsize=11, color="#888")

fig.suptitle(
    "Wet-Scavenging & Rebound Analysis — Chandrapur Monsoon 2026\n"
    "Events detected from coherent network-wide PM₂.₅ drops (≥5 stations, >20% fall in 3h)",
    fontsize=12.5, fontweight="bold", y=1.01
)

plt.savefig(FIG_OUT, dpi=200, bbox_inches="tight")
print(f"\nSaved → {FIG_OUT}")
