import csv

col_totals = None
n_rows = 0

with open('power_scan.csv') as f:
    for row in csv.reader(f):
        row = [x.strip() for x in row]
        vals = [float(x) for x in row[6:]]
        if col_totals is None:
            col_totals = [0.0] * len(vals)
        for i, v in enumerate(vals):
            col_totals[i] += v
        n_rows += 1

avg = [t / n_rows for t in col_totals]
peak_idx = avg.index(max(avg))
low_hz = 54280000
high_hz = 54340000
step_hz = (high_hz - low_hz) / len(avg)
peak_freq = low_hz + peak_idx * step_hz
overall_avg = sum(avg) / len(avg)

print(f"Rows processed: {n_rows}")
print(f"Peak averaged bin: index {peak_idx}, ~{peak_freq/1e6:.4f} MHz, avg power {avg[peak_idx]:.2f} dB")
print(f"Overall average power: {overall_avg:.2f} dB")
print(f"Peak stands out by: {avg[peak_idx] - overall_avg:.2f} dB above average")
