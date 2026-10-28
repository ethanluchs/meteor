import pandas as pd

df = pd.read_csv('events_overnight.csv')
print(f"Total events: {len(df)}\n")

if len(df) > 0:
    df['start_iso'] = pd.to_datetime(df['start_iso'])
    df['hour'] = df['start_iso'].dt.hour
    counts = df.groupby('hour').size().reindex(range(24), fill_value=0)
    print("=== Events by hour ===")
    print(counts)
    print("\n=== Full event list ===")
    print(df[['start_iso', 'duration_s', 'peak_snr', 'dominant_bin', 'dominant_offset_hz', 'dominant_fraction']].to_string())
