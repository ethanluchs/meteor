import pandas as pd

df = pd.read_csv('raw_snr_log.csv')
df['dt'] = pd.to_datetime(df['unix_time'], unit='s')

# per-row excess above the tracked baseline, in units of its own std
df['excess_sigma'] = (df['snr_db'] - df['baseline_mean']) / df['baseline_std']

threshold = 2.0  # looser than a single-frame test needs to be, since we're
                 # about to require sustained elevation across many rows
df['elevated'] = df['excess_sigma'] > threshold

# group consecutive elevated rows into runs, tolerating small gaps
gap_tolerance = 3  # rows (~0.5s at this sampling rate)
df['gap'] = (~df['elevated']).astype(int)
df['gap_run'] = df['gap'].groupby((df['elevated']).cumsum()).cumsum()
df['bridge'] = df['elevated'] | ((df['gap_run'] <= gap_tolerance) & df['elevated'].shift(1).fillna(False))

df['group'] = (df['elevated'] != df['elevated'].shift()).cumsum()
runs = df[df['elevated']].groupby('group').agg(
    start=('dt', 'first'),
    end=('dt', 'last'),
    duration_s=('dt', lambda x: (x.max() - x.min()).total_seconds()),
    n_rows=('dt', 'count'),
    peak_snr=('snr_db', 'max'),
    dominant_bin=('peak_bin', lambda x: x.mode()[0]),
).reset_index(drop=True)

real_events = runs[runs['duration_s'] >= 1.0].sort_values('start')
print(f"Total elevated runs: {len(runs)}")
print(f"Runs lasting 1+ second: {len(real_events)}\n")
print(real_events.to_string())

if len(real_events) > 0:
    real_events['hour'] = real_events['start'].dt.hour
    print("\n=== Events by hour ===")
    print(real_events.groupby('hour').size().reindex(range(24), fill_value=0))
