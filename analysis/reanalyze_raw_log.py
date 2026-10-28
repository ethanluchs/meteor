#!/usr/bin/env python3
"""
Re-analyze raw_snr_log.csv using a windowed bin-occupancy statistic
(the same idea as the old mode_fraction) instead of a per-frame SNR threshold.

Key differences from the previous analysis:
  1. Time is reconstructed from the row counter and the known frame rate,
     not from the 1-second-resolution unix_time column.
  2. The detection statistic is "what fraction of frames in this window had
     their peak at the carrier bin", which averages down as 1/sqrt(N).
     Per-frame snr_db does not.
  3. The significance bar is set empirically from decoy bins in the same
     window, so it does not depend on any Gaussian/independence assumption.

Usage:  python3 reanalyze_raw_log.py [path/to/raw_snr_log.csv]
"""

import sys
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- config ----
CSV            = sys.argv[1] if len(sys.argv) > 1 else 'raw_snr_log.csv'
SAMPLE_RATE    = 250000.0   # Hz
FFT_SIZE       = 4096
LOG_DECIMATION = 10         # raw log writes one row per Nth frame
BIN_TOL        = 1          # carrier +/- this many bins counts as a "hit"
WINDOW_S       = 10.0       # detection window length, seconds
N_DECOYS       = 32         # control bins used to build the empirical null

FRAME_RATE = SAMPLE_RATE / FFT_SIZE          # ~61.035 frames/s
ROW_DT     = LOG_DECIMATION / FRAME_RATE     # ~0.16384 s between logged rows
WIN        = max(3, int(round(WINDOW_S / ROW_DT)))

# ------------------------------------------------------------------ load ----
df = pd.read_csv(CSV)
n = len(df)
if n < 10 * WIN:
    sys.exit(f"only {n} rows, need at least {10 * WIN}")

t0 = int(df['unix_time'].iloc[0])
t_rel = np.arange(n, dtype=float) * ROW_DT
df['dt'] = pd.to_datetime(t0, unit='s') + pd.to_timedelta(t_rel, unit='s')

clock_span   = int(df['unix_time'].iloc[-1]) - t0
counter_span = t_rel[-1]
print("=" * 72)
print("TIMEBASE")
print("=" * 72)
print(f"rows                 : {n}")
print(f"row spacing          : {ROW_DT * 1000:.2f} ms  ({1/ROW_DT:.2f} rows/s)")
print(f"span from row counter: {counter_span / 3600:.4f} h")
print(f"span from unix_time  : {clock_span / 3600:.4f} h")
print(f"drift                : {clock_span - counter_span:+.1f} s "
      f"({100 * (clock_span - counter_span) / max(counter_span, 1):+.3f} %)")
print("  (a large positive drift means frames were dropped / the loop stalled;")
print("   a large negative drift means FRAME_RATE or LOG_DECIMATION is wrong)")

# ----------------------------------------------------- find carrier bin ----
counts = df['peak_bin'].value_counts()
lo, hi = int(df['peak_bin'].min()), int(df['peak_bin'].max())
n_bins = hi - lo + 1
carrier = int(counts.index[0])

print()
print("=" * 72)
print("BIN OCCUPANCY (whole file)")
print("=" * 72)
print(f"search window spans bins {lo}..{hi}  ({n_bins} bins, "
      f"{n_bins * SAMPLE_RATE / FFT_SIZE / 1000:.1f} kHz)")
print(f"uniform-noise share per bin: {1 / n_bins:.4f}")
print("\ntop 12 bins:")
for b, c in counts.head(12).items():
    print(f"  bin {int(b):4d}  offset {(int(b) - carrier) * SAMPLE_RATE / FFT_SIZE:+8.1f} Hz"
          f"  {c:7d} rows  {c / n:.4f}  ({c / n * n_bins:5.1f}x uniform)")
print(f"\n-> carrier bin taken as {carrier}")


def hit_series(center):
    return ((df['peak_bin'] - center).abs() <= BIN_TOL).to_numpy(dtype=float)


hit = hit_series(carrier)
p0 = hit.mean()
uniform = (2 * BIN_TOL + 1) / n_bins
print(f"   carrier hit rate (+/-{BIN_TOL} bins): {p0:.4f}")
print(f"   uniform-noise expectation          : {uniform:.4f}")
print(f"   ratio                              : {p0 / uniform:.1f}x")
if p0 / uniform < 2:
    print("   WARNING: carrier is not measurably winning more than a random bin.")
    print("   Either the tuning/offset is wrong or the carrier is below the floor.")

# ----------------------------------------------- windowed statistic + null --
def windowed_z(center):
    """Rolling count of hits on `center`, converted to a binomial z-score."""
    h = pd.Series(hit_series(center))
    k = h.rolling(WIN, center=True, min_periods=WIN).sum()
    p = h.mean()
    if p <= 0:
        return k, np.full(n, np.nan), p
    mu = WIN * p
    sd = np.sqrt(WIN * p * (1 - p))
    return k.to_numpy(), ((k - mu) / sd).to_numpy(), p


k_car, z_car, _ = windowed_z(carrier)

rng = np.random.default_rng(0)
avoid = set(range(carrier - 4 * (BIN_TOL + 1), carrier + 4 * (BIN_TOL + 1) + 1))
pool = [b for b in range(lo + BIN_TOL, hi - BIN_TOL + 1) if b not in avoid]
decoys = rng.choice(pool, size=min(N_DECOYS, len(pool)), replace=False)

decoy_max = []
for b in decoys:
    _, z, p = windowed_z(int(b))
    if np.isfinite(z).any():
        decoy_max.append(np.nanmax(z))
decoy_max = np.array(decoy_max)

z_bar = float(np.nanmax(decoy_max)) if len(decoy_max) else 5.0

print()
print("=" * 72)
print(f"WINDOWED DETECTION  (window = {WIN} rows = {WIN * ROW_DT:.1f} s)")
print("=" * 72)
print(f"carrier: peak windowed z = {np.nanmax(z_car):.2f}  "
      f"(max hits in a window: {int(np.nanmax(k_car))}/{WIN} = "
      f"{np.nanmax(k_car) / WIN:.3f})")
print(f"decoys : {len(decoy_max)} control bins, peak z per bin "
      f"median {np.median(decoy_max):.2f}, max {z_bar:.2f}")
print(f"\nbar for a real detection: carrier z must beat every decoy's best "
      f"window (z > {z_bar:.2f})")

if np.nanmax(z_car) <= z_bar:
    print("\nRESULT: nothing in this run beats the empirical noise floor.")
    print("The night contains no windowed carrier enhancement that a random")
    print("bin in the same data does not also produce. This is a clean null,")
    print("not a failed analysis.")
else:
    print("\nRESULT: at least one window beats every decoy bin. Candidates below.")

# ------------------------------------------------------- list candidates ----
cand = pd.Series(z_car > z_bar).fillna(False).to_numpy()
if cand.any():
    grp = np.cumsum(np.concatenate(([0], np.diff(cand.astype(int)) != 0)))
    out = []
    for g in np.unique(grp[cand]):
        idx = np.where(cand & (grp == g))[0]
        a, b = idx[0] - WIN // 2, idx[-1] + WIN // 2
        a, b = max(a, 0), min(b, n - 1)
        out.append({
            'start': df['dt'].iloc[a],
            'duration_s': round((b - a) * ROW_DT, 2),
            'peak_z': round(float(np.nanmax(z_car[idx])), 2),
            'peak_frac': round(float(np.nanmax(k_car[idx]) / WIN), 3),
            'base_frac': round(p0, 3),
            'mean_snr': round(float(df['snr_db'].iloc[a:b + 1].mean()), 2),
            'max_snr': round(float(df['snr_db'].iloc[a:b + 1].max()), 2),
        })
    ev = pd.DataFrame(out).sort_values('peak_z', ascending=False)
    print()
    print(ev.to_string(index=False))

    ev['hour'] = pd.to_datetime(ev['start']).dt.hour
    hours = sorted(df['dt'].dt.hour.unique())
    print("\nevents per hour (only hours actually covered by this run):")
    print(ev.groupby('hour').size().reindex(hours, fill_value=0).to_string())

# ------------------------------------------------- baseline sanity check ----
print()
print("=" * 72)
print("BASELINE TRACKER HEALTH")
print("=" * 72)
bm = df['baseline_mean']
print(f"baseline_mean: {bm.min():.2f} .. {bm.max():.2f} dB, "
      f"std over the run {bm.std():.3f} dB")
print(f"baseline_std : mean {df['baseline_std'].mean():.3f} dB")
excess = (df['snr_db'] - bm) / df['baseline_std']
print(f"per-frame excess sigma: mean {excess.mean():+.3f}, "
      f"p99 {excess.quantile(0.99):.2f}, max {excess.max():.2f}")
print(f"fraction above 2 sigma: {(excess > 2).mean():.4f}  "
      f"(Gaussian would give 0.0228)")
print(f"fraction above 3 sigma: {(excess > 3).mean():.6f}  "
      f"(Gaussian would give 0.001350)")
print(f"fraction above 4 sigma: {(excess > 4).mean():.6f}  "
      f"(Gaussian would give 0.000032)")
print("\nIf the 3/4-sigma fractions are far BELOW Gaussian, the EMA baseline is")
print("chasing the signal and eating real excursions. If they are far ABOVE,")
print("the per-frame SNR distribution is right-skewed (expected for a")
print("max-over-bins statistic) and sigma thresholds are not calibrated.")
