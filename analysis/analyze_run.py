import pandas as pd

df = pd.read_csv('carrier_data.csv')
df['iso_time'] = pd.to_datetime(df['iso_time'])

print(f"Total rows: {len(df)}")
print(f"Time range: {df['iso_time'].min()} to {df['iso_time'].max()}\n")

print("=== Overall stats ===")
print(df[['noise_db', 'snr_median', 'snr_max', 'mode_fraction']].describe())

print("\n=== Most common mode_bin (top 10) ===")
print(df['mode_bin'].value_counts().head(10))

print("\n=== Rows with unusually high mode_fraction (top 20) ===")
print(df.nlargest(20, 'mode_fraction')[['iso_time', 'mode_bin', 'mode_offset_hz', 'mode_fraction', 'snr_max']])

print("\n=== Rows with unusually high snr_max (top 20) ===")
print(df.nlargest(20, 'snr_max')[['iso_time', 'mode_bin', 'mode_offset_hz', 'snr_max', 'mode_fraction']])

# Flag candidate "events": mode_fraction above a threshold, meaning one bin
# dominated an unusually large share of frames in that bucket
threshold = 0.10
candidates = df[df['mode_fraction'] > threshold].sort_values('iso_time')
print(f"\n=== Candidate events (mode_fraction > {threshold}): {len(candidates)} rows ===")
print(candidates[['iso_time', 'mode_bin', 'mode_offset_hz', 'mode_fraction', 'snr_max']].to_string())
