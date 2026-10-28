import pandas as pd

df = pd.read_csv('carrier_data.csv')
df['iso_time'] = pd.to_datetime(df['iso_time'])

# Exclude the suspected cluster window when computing the "quiet" baseline,
# so the event doesn't inflate its own threshold
quiet = df[~((df['iso_time'] >= '2026-09-18 20:37:00') & (df['iso_time'] <= '2026-09-18 21:17:00'))]

mean_mf = quiet['mode_fraction'].mean()
std_mf = quiet['mode_fraction'].std()
threshold = mean_mf + 5 * std_mf

print(f"Baseline mean: {mean_mf:.4f}, std: {std_mf:.4f}")
print(f"5-sigma threshold: {threshold:.4f}")

candidates = df[df['mode_fraction'] > threshold].sort_values('iso_time')
print(f"\nCandidate events above threshold: {len(candidates)}")
print(candidates[['iso_time', 'mode_bin', 'mode_offset_hz', 'mode_fraction', 'snr_max']].to_string())
