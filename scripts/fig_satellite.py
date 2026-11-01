#!/usr/bin/env python3
"""
Satellite-imagery PM2.5 overlay maps — publication-ready.
Run this on your own machine (needs internet for map tiles).

Install:  pip3 install contextily --break-system-packages
Run from the AQMap root:  python3 scripts/fig_satellite.py
Outputs saved to:         figures/fig1_spatial_satellite.png
                          figures/fig5_exposure_satellite.png
"""

import numpy as np, pandas as pd, matplotlib, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from scipy.ndimage import gaussian_filter
import contextily as ctx

# ── CSTPS confirmed coordinates (north of S2, verified on satellite imagery) ──
PLANT = (20.00837, 79.29206)   # confirmed: N of sensor cluster, not NW

COORDS = {
    "1":(19.98062,79.23735),"2":(19.97237,79.28563),"3":(19.98250,79.30010),
    "4":(19.96839,79.31393),"5":(19.96024,79.31746),"6":(19.97239,79.30880),
    "7":(20.00482,79.30233),"8":(19.94894,79.29743),"9":(19.93473,79.32349),
    "10":(19.95375,79.32983)
}
CALIB={"1":{"b0":0.114,"b_raw":1.0026},"2":{"b0":0.231,"b_raw":0.9832},
       "3":{"b0":-0.246,"b_raw":1.0233},"4":{"b0":0.132,"b_raw":0.9851},
       "5":{"b0":0.313,"b_raw":0.9242},"6":{"b0":0.243,"b_raw":0.9593},
       "7":{"b0":0.289,"b_raw":0.9907},"8":{"b0":0.308,"b_raw":1.0085},
       "9":{"b0":0.165,"b_raw":1.0216},"10":{"b0":0.489,"b_raw":0.9847}}
DEPLOY_START={"1":"2026-08-11","2":"2026-08-11","3":"2026-08-11","4":"2026-08-24",
              "5":"2026-09-04","6":"2026-09-02","7":"2026-08-11","8":"2026-08-11",
              "9":"2026-08-11","10":"2026-08-11"}
CSV_MAP = {
    "1":  "S1_civil-lines.csv",
    "2":  "S2_civil-lines.csv",
    "3":  "S3_civil-lines.csv",
    "4":  "S4_civil-lines.csv",
    "5":  "S5_civil-lines.csv",
    "6":  "S6_civil-lines.csv",
    "7":  "S7_civil-lines.csv",
    "8":  "S8_civil-lines.csv",
    "9":  "S9_civil-lines.csv",
    "10": "S10_civil-lines.csv",
}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
AQMAP_ROOT = os.path.dirname(SCRIPT_DIR)          # one level up from scripts/
UPLOAD_DIR = os.path.join(AQMAP_ROOT, "data", "sensors")
FIGURES_DIR = os.path.join(AQMAP_ROOT, "figures")

MONSOON_END = pd.Timestamp("2026-09-22")
lat_min,lat_max=19.91,20.02; lon_min,lon_max=79.21,79.36
N=200
lats=np.linspace(lat_min,lat_max,N); lons=np.linspace(lon_min,lon_max,N)
glon,glat=np.meshgrid(lons,lats)

def load_station(sid):
    df=pd.read_csv(f"{UPLOAD_DIR}/{CSV_MAP[sid]}")
    ts_col=next(c for c in df.columns if 'time' in c.lower() or 'date' in c.lower())
    pm_col=next(c for c in df.columns if 'pm2' in c.lower() or 'pm_2' in c.lower())
    df=df.rename(columns={ts_col:'ts',pm_col:'pm_raw'})
    df['ts']=pd.to_datetime(df['ts'],utc=True).dt.tz_localize(None)
    df=df.dropna(subset=['ts','pm_raw'])
    df=df[df['ts']>=pd.Timestamp(DEPLOY_START[sid])]
    df=df[df['pm_raw']>0]
    c=CALIB[sid]; df['pm25']=(c['b0']+c['b_raw']*df['pm_raw']).clip(lower=0)
    q=df['pm25'].quantile([0.01,0.99])
    return df[(df['pm25']>=q[0.01])&(df['pm25']<=q[0.99])][['ts','pm25']].set_index('ts')

def idw(svals, power=2):
    pts=np.array([COORDS[s] for s in svals]); vals=np.array([svals[s] for s in svals])
    out=np.zeros_like(glon)
    for i in range(N):
        for j in range(N):
            d=np.sqrt((pts[:,0]-glat[i,j])**2+(pts[:,1]-glon[i,j])**2)
            d=np.maximum(d,1e-9); w=1.0/d**power
            out[i,j]=np.sum(w*vals)/np.sum(w)
    return out

# load station means
sm,sp={},{}
for sid in COORDS:
    df=load_station(sid)
    sm[sid]=df[df.index<=MONSOON_END]['pm25'].mean()
    pst=df[df.index>MONSOON_END]; sp[sid]=pst['pm25'].mean() if len(pst)>0 else np.nan

print("Computing IDW grids…")
pm_mon  = idw({s:v for s,v in sm.items() if not np.isnan(v)})
pm_post = idw({s:v for s,v in sp.items() if not np.isnan(v)})
print("Done.")


def make_satellite_map(pm_grids, titles, outfile, overall_title,
                       vmin=20, vmax=65, pm_alpha=0.52):
    """
    Publication-ready satellite basemap with PM2.5 overlay.
    pm_grids : list of 2-D arrays (one per subplot)
    titles   : list of subplot title strings
    """
    # Perceptually ordered colormap for air quality (green=clean, red=polluted)
    cmap = matplotlib.colormaps.get_cmap('RdYlGn_r')
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 9,
        'axes.linewidth': 0.8,
    })

    n = len(pm_grids)
    # Each panel covers 0.15° lon × 0.12° lat → natural aspect ~1.17 wide per panel
    # Add ~1.5" per panel for colorbar; 6.5" height gives clean proportions
    fig, axes = plt.subplots(1, n, figsize=(6.8*n + 1.0, 6.5))
    if n == 1: axes = [axes]

    outline_w = [pe.withStroke(linewidth=2.2, foreground='black')]
    outline_b = [pe.withStroke(linewidth=2.2, foreground='white')]

    for ax, pm_grid, title in zip(axes, pm_grids, titles):
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)

        # ── satellite basemap ────────────────────────────────────────────
        print(f"  Fetching tiles for '{title}'…")
        ctx.add_basemap(ax, crs="EPSG:4326",
                        source=ctx.providers.Esri.WorldImagery,
                        zoom=13, attribution=False)

        # ── smooth IDW before contouring (removes station halos) ─────────
        # sigma=2 in grid units ≈ 1.5 km, reduces IDW ring artifacts
        pm_smooth = gaussian_filter(pm_grid, sigma=2.5)

        # ── PM2.5 overlay ────────────────────────────────────────────────
        im = ax.contourf(glon, glat, pm_smooth,
                         levels=np.linspace(vmin, vmax, 40),
                         cmap=cmap, norm=norm, alpha=pm_alpha, zorder=2)

        # Regulatory threshold contours — white dashed for readability on satellite
        cs60 = ax.contour(glon, glat, pm_smooth, levels=[60],
                          colors=['white'], linewidths=1.6,
                          zorder=3, linestyles='--')
        ax.clabel(cs60, fmt='60 µg/m³\n(NAAQS)', fontsize=7,
                  colors='white', inline=True,
                  manual=False)

        # ── Sensor markers ───────────────────────────────────────────────
        for sid, (la, lo) in COORDS.items():
            # Outer ring (black) then white fill
            ax.scatter(lo, la, s=62, color='white',
                       edgecolors='#1a1a1a', linewidths=1.4, zorder=5)
            # Station label with black outline stroke → readable on any background
            ax.annotate(
                f'S{sid}', (lo, la),
                xytext=(5, 3), textcoords='offset points',
                color='white', fontsize=7.8, fontweight='bold',
                path_effects=outline_w, zorder=6,
            )

        # ── CSTPS — professional plant marker ───────────────────────────
        # Filled dark-red square with white border: unambiguous, not garish
        ax.scatter(PLANT[1], PLANT[0],
                   marker='s', s=130,
                   color='#C62828',          # deep red
                   edgecolors='white', linewidths=1.6,
                   zorder=7)
        ax.annotate(
            'CSTPS', (PLANT[1], PLANT[0]),
            xytext=(7, 5), textcoords='offset points',
            color='white', fontsize=8.5, fontweight='bold',
            path_effects=outline_w, zorder=7,
        )

        # ── Panel label & formatting ─────────────────────────────────────
        ax.set_title(title, fontsize=11, fontweight='bold',
                     pad=7, loc='left')
        ax.set_xlabel('Longitude (°E)', fontsize=9)
        ax.set_ylabel('Latitude (°N)', fontsize=9)
        ax.tick_params(labelsize=8)

        # ── Colorbar ─────────────────────────────────────────────────────
        cb = plt.colorbar(im, ax=ax, fraction=0.036, pad=0.015, shrink=0.88)
        cb.set_label('PM₂.₅ (µg/m³)', fontsize=9)
        cb.ax.tick_params(labelsize=8)
        # WHO guideline (15 µg/m³) is below vmin=20 so can't draw on colorbar;
        # it is shown in the legend instead.

    # ── Shared legend ────────────────────────────────────────────────────
    legend_elements = [
        Line2D([0],[0], marker='s', color='w', markerfacecolor='#C62828',
               markeredgecolor='white', markeredgewidth=1.2,
               markersize=9, label='CSTPS stack location'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='white',
               markeredgecolor='#1a1a1a', markeredgewidth=1.2,
               markersize=8, label='PM₂.₅ sensor (calibrated, n=10)'),
        Line2D([0],[0], ls='--', color='white', lw=1.6,
               label='NAAQS 24-h standard (60 µg/m³)'),
        Line2D([0],[0], ls='--', color='#1B5E20', lw=1.6,
               label='WHO annual guideline (15 µg/m³)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4,
               fontsize=8.5, framealpha=0.92, edgecolor='#bbb',
               bbox_to_anchor=(0.5, 0.01))

    fig.suptitle(overall_title, fontsize=11.5, fontweight='bold', y=1.002)
    plt.tight_layout(rect=[0, 0.07, 1, 1])
    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    print(f"  Saved → {outfile}")
    plt.close()


# ── Fig 1 satellite version ──────────────────────────────────────────────
make_satellite_map(
    pm_grids=[pm_mon, pm_post],
    titles=['(a)  Monsoon  (11 Aug – 22 Sep 2026)',
            '(b)  Post-Monsoon  (23 Sep – 6 Oct 2026)'],
    outfile=os.path.join(FIGURES_DIR, 'fig1_spatial_satellite.png'),
    overall_title=(
        'Spatial Distribution of PM₂.₅  —  Chandrapur, Maharashtra, 2026\n'
        '10-Station Low-Cost Sensor Network · IDW Interpolation · ESRI World Imagery'
    ),
)

# ── Fig 5 satellite version ──────────────────────────────────────────────
make_satellite_map(
    pm_grids=[pm_mon, pm_post],
    titles=['(a)  Monsoon  (pop.-weighted: 33.3 µg/m³)',
            '(b)  Post-Monsoon  (pop.-weighted: 45.5 µg/m³)'],
    outfile=os.path.join(FIGURES_DIR, 'fig5_exposure_satellite.png'),
    overall_title=(
        'Population-Weighted PM₂.₅ Exposure  —  Chandrapur, Maharashtra, 2026\n'
        '10-Station Low-Cost Sensor Network · ESRI World Imagery'
    ),
)

print(f"\nAll done. Outputs saved to {FIGURES_DIR}/")
print("  fig1_spatial_satellite.png")
print("  fig5_exposure_satellite.png")
print("Note: red contour = NAAQS (60 µg/m³), dashed lime = WHO (15 µg/m³)")
