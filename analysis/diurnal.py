import pandas as pd

df = pd.read_csv('events.csv')
df['start_iso'] = pd.to_datetime(df['start_iso'])
df['hour'] = df['start_iso'].dt.hour

counts = df.groupby('hour').size().reindex(range(24), fill_value=0)
print(counts)
