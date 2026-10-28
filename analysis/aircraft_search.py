import pandas as pd

df = pd.read_csv('carrier_run2.csv')
df['iso_time'] = pd.to_datetime(df['iso_time'])

mean_mf = df['mode_fraction'].mean()
std_mf = df['mode_fraction'].std()
threshold = mean_mf + 3 * std_mf  # a bit looser than 5-sigma, to catch longer/weaker events too

df['elevated'] = df['mode_fraction'] > threshold
df['group'] = (df['elevated'] != df['elevated'].shift()).cumsum()

runs = df[df['elevated']].groupby('group').agg(
    start=('iso_time', 'first'),
    end=('iso_time', 'last'),
    duration_s=('iso_time', lambda x: (x.max() - x.min()).total_seconds() + 1),
    n_rows=('iso_time', 'count'),
    offset_start=('mode_offset_hz', 'first'),
    offset_end=('mode_offset_hz', 'last'),
).reset_index(drop=True)

# Aircraft passes tend to run noticeably longer than meteor pings
long_runs = runs[runs['duration_s'] >= 20].sort_values('duration_s', ascending=False)
print(f"Total elevated runs found: {len(runs)}")
print(f"Runs lasting 20+ seconds: {len(long_runs)}\n")
print(long_runs.to_string())

print("\n=== Duration distribution across ALL runs ===")
print(runs['duration_s'].describe())
print("\n=== Top 15 longest runs regardless of cutoff ===")
print(runs.sort_values('duration_s', ascending=False).head(15).to_string())
