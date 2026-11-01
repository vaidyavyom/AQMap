#!/usr/bin/env python3
"""
Reference-anchored calibration check + fit  (Chauhan Colony CPCB ref, node R1).

Run this AFTER a couple of days of the poller logging both your network and R1 into
history.csv.  It:

  1. time-matches your network (hourly median) against the reference,
  2. tests whether your excess over the reference grows with humidity
     (the signature of hygroscopic over-reading),
  3. proposes a correction that ties your network to the reference BUT is biased so
     it rarely reads *below* the reference (your "don't underestimate" rule).

It only PRINTS suggested coefficients - it does not touch config.py. Review the
numbers, then paste the CALIB block in if you're happy.

    python3 ref_calibration.py
"""
import csv, statistics as st, datetime as dt
from config import CALIB

HISTORY   = "history.csv"
REF_ID    = "R1"
MIN_HOURS = 24        # need at least this many matched hours to say anything
LAGS      = [0, 1, 2] # test aligning network hour t with reference hour t-lag (ref lags ~1-2h)
NO_UNDERESTIMATE_PCTL = 10   # keep >=90% of corrected values at/above the reference

# ---------- tiny stats helpers ----------
def ols(xs, ys):
    """y ~ b0 + b1*x ; returns (b1, b0, r)."""
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); syy=sum((y-my)**2 for y in ys)
    sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    if sxx==0 or syy==0: return 0.0, my, 0.0
    b1=sxy/sxx; return b1, my-b1*mx, sxy/(sxx*syy)**0.5
def pctl(v, p):
    s=sorted(v); k=(len(s)-1)*p/100; f=int(k)
    return s[f] if f+1>=len(s) else s[f]+(s[f+1]-s[f])*(k-f)

# ---------- load history, aggregate to hourly ----------
def load():
    ref={}     # hour -> list of reference pm25
    net={}     # hour -> list of (unit calibrated pm25)
    netrh={}   # hour -> list of unit RH
    for r in csv.DictReader(open(HISTORY)):
        try: h=dt.datetime.fromisoformat(r["poll_time"]).replace(minute=0,second=0,microsecond=0)
        except Exception: continue
        mon=r["monitor"]
        try: cal=float(r["calibrated_pm25"])
        except Exception: continue
        if mon==REF_ID:
            ref.setdefault(h,[]).append(cal)
        else:
            net.setdefault(h,[]).append(cal)
            try: netrh.setdefault(h,[]).append(float(r["humidity"]))
            except Exception: pass
    refh={h:st.mean(v) for h,v in ref.items()}
    neth={h:st.median(v) for h,v in net.items()}
    rhh ={h:st.median(v) for h,v in netrh.items() if v}
    return refh, neth, rhh

def pair(refh, neth, rhh, lag):
    """align network hour t with reference hour (t - lag)."""
    rows=[]
    for h,nv in neth.items():
        rh=rhh.get(h)
        rk=h - dt.timedelta(hours=lag)
        if rk in refh and rh is not None:
            rows.append((h, nv, refh[rk], rh))
    return rows

def main():
    refh, neth, rhh = load()
    if not refh:
        print("No reference (R1) rows in history.csv yet. Let the poller run, then re-run."); return
    # pick the lag that best correlates net vs ref
    best=None
    for lag in LAGS:
        rows=pair(refh,neth,rhh,lag)
        if len(rows)<6: continue
        _,_,r=ols([x[1] for x in rows],[x[2] for x in rows])
        if best is None or r>best[0]: best=(r,lag,rows)
    if best is None:
        print("Not enough overlapping hours yet. Re-run after more data."); return
    r_best, lag, rows = best
    n=len(rows)
    net=[x[1] for x in rows]; ref=[x[2] for x in rows]; rh=[x[3] for x in rows]
    bias=st.mean([a-b for a,b in zip(net,ref)])
    print(f"matched hours: {n}   (best ref lag = {lag} h, r = {r_best:.2f})")
    print(f"network median mean: {st.mean(net):.1f}   reference mean: {st.mean(ref):.1f}   "
          f"mean excess (net-ref): {bias:+.1f} ug/m3")
    if n < MIN_HOURS:
        print(f"\n--> only {n} hours; need >= {MIN_HOURS} for a trustworthy fit. Re-run later.")
        return

    # does the excess grow with humidity?
    excess=[a-b for a,b in zip(net,ref)]
    b_rh, a_rh, r_rh = ols(rh, excess)
    print(f"\nexcess vs RH:  excess ~= {a_rh:+.2f} {b_rh:+.3f}*RH   (r = {r_rh:.2f})")
    if b_rh>0 and r_rh>0.3:
        print("  -> excess RISES with humidity: consistent with hygroscopic over-reading.")
    else:
        print("  -> no clear humidity trend: the gap is mostly spatial, not an RH artifact.")

    # centered correction: final = net - (a_rh + b_rh*RH)  ties net onto ref on average
    centered=[nv-(a_rh+b_rh*rv) for nv,rv in zip(net,rh)]
    resid=[c-rf for c,rf in zip(centered,ref)]          # want most of these >= 0
    margin=max(0.0, -pctl(resid, NO_UNDERESTIMATE_PCTL))  # lift so ~90% stay >= ref
    print(f"\nconservative margin (+{margin:.1f} ug/m3) keeps ~{100-NO_UNDERESTIMATE_PCTL}% of readings "
          f"at or above the reference.")

    # express as per-unit CALIB update: final = (b0 - a_rh + margin) + b_raw*raw + (-b_rh)*RH
    print("\nSuggested CALIB (reference-anchored, no-underestimate). Review, then paste into config.py:\n")
    print("CALIB = {")
    for uid in sorted(CALIB, key=lambda x:(len(x),x)):
        c=CALIB[uid]; b0=c.get("b0",0.0); braw=c.get("b_raw",1.0)
        nb0=round(b0 - a_rh + margin, 3); nbrh=round(-b_rh, 4)
        print(f'    "{uid}": {{"b0": {nb0}, "b_raw": {braw}, "b_rh": {nbrh}}},')
    print("}")
    print("\n(To also apply this to the map, keep RH_CORRECTION as-is; the b_rh term above "
          "is a linear RH correction inside apply_calibration, separate from the kappa knob.)")

if __name__ == "__main__":
    main()
