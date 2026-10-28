import numpy as np, pandas as pd

BIN_HZ = 250000.0/16384; DT = 16384/250000.0
df = pd.read_csv('raw_snr_log.csv')
df['locked'] = (df.peak_bin - df.carrier_bin).abs() <= 3
grp = (df.locked != df.locked.shift()).cumsum()

print(f"rows {len(df)}, carrier-locked {df.locked.mean():.3f}")
print(f"snr_db: locked median {df[df.locked].snr_db.median():.2f}, "
      f"unlocked median {df[~df.locked].snr_db.median():.2f}\n")

runs = []
for g, d in df[df.locked].groupby(grp):
    if len(d) < 15: continue
    t = d.t_sec.values - d.t_sec.values[0]
    y = d.frac_bin.values
    if t[-1] <= 0: continue
    slope, icept = np.polyfit(t, y, 1)
    resid = y - (slope*t + icept)
    ss = 1 - resid.var()/max(y.var(), 1e-12)
    s = d.snr_db.values
    lo, hi = s.min(), s.max()
    thr10, thr90 = lo + .1*(hi-lo), lo + .9*(hi-lo)
    try:
        rise = t[np.argmax(s >= thr90)] - t[np.argmax(s >= thr10)]
    except Exception:
        rise = np.nan
    runs.append(dict(t0=d.t_sec.values[0], dur=t[-1], n=len(d),
                     drift_hz=slope*BIN_HZ*t[-1], hz_per_s=slope*BIN_HZ,
                     r2=ss, snr_max=hi, snr_rise=hi-lo, rise_s=rise))
r = pd.DataFrame(runs)
print(f"carrier-locked runs >=1s: {len(r)}\n")

print("=== AIRCRAFT TEST: smooth monotonic drift ===")
a = r[(r.dur >= 3) & (r.r2 > 0.6) & (r.drift_hz.abs() > 8)]
print(a.sort_values('drift_hz', key=abs, ascending=False).head(10).to_string(index=False)
      if len(a) else "none — no smooth frequency sweeps found")

print("\n=== METEOR TEST: fast rise ===")
m = r[(r.snr_rise > 4) & (r.rise_s < 1.0)]
print(m.sort_values('snr_rise', ascending=False).head(10).to_string(index=False)
      if len(m) else "none — nothing rises in under a second")

print("\n=== slowest-to-fastest rise, strongest 10 runs ===")
print(r.nlargest(10, 'snr_rise')[['t0','dur','snr_max','snr_rise','rise_s','drift_hz','r2']].to_string(index=False))

print("\n=== baseline drift over the run (10-min blocks) ===")
df['blk'] = (df.t_sec//600).astype(int)
print(df.groupby('blk').agg(locked_frac=('locked','mean'),
                            snr_med=('snr_db','median')).to_string())
