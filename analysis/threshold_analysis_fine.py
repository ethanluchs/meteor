import pandas as pd

df = pd.read_csv('carrier_run2.csv')
df['iso_time'] = pd.to_datetime(df['iso_time'])

print(f"Total rows: {len(df)}")
print(f"Time range: {df['iso_time'].min()} to {df['iso_time'].max()}\n")

print("=== Overall stats ===")
print(df[['noise_db', 'snr_median', 'snr_max', 'mode_fraction']].describe())

print("\n=== Most common mode_bin (top 10) ===")
print(df['mode_bin'].value_counts().head(10))

mean_mf = df['mode_fraction'].mean()
std_mf = df['mode_fraction'].std()
threshold = mean_mf + 5 * std_mf
print(f"\nBaseline mean: {mean_mf:.4f}, std: {std_mf:.4f}")
print(f"5-sigma threshold: {threshold:.4f}")

candidates = df[df['mode_fraction'] > threshold].sort_values('iso_time')
print(f"\nCandidate events above threshold: {len(candidates)}")
print(candidates[['iso_time', 'mode_bin', 'mode_offset_hz', 'mode_fraction', 'snr_max']].head(30).to_string())
