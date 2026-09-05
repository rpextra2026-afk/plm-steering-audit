"""
Quick free check enabled by the newly-uploaded CSV: is collapse in the orthogonalized-foldability
run confounded by sequence length? Uses analysis/orthogonalized_foldability_and_refold_sequences.csv
(from 39-ai4dd-orthogonalized-foldability-and-refold-fix.ipynb). No GPU, no new experiment.
"""
import csv
from collections import defaultdict

path = "orthogonalized_foldability_and_refold_sequences.csv"

by_cond = defaultdict(list)
with open(path, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        by_cond[row["condition"]].append(row)

print(f"{'Condition':22s} {'N':>4s} {'mean len':>9s} {'mean usable':>12s} {'mean pLDDT':>11s} {'collapse%':>10s}")
print("-" * 74)
for cond, rows in by_cond.items():
    n = len(rows)
    mean_len = sum(float(r["gen_length"]) for r in rows) / n
    mean_usable = sum(float(r["usable_length"]) for r in rows) / n
    mean_plddt = sum(float(r["plddt"]) for r in rows) / n
    collapse = sum(float(r["collapse"]) for r in rows) / n
    print(f"{cond:22s} {n:4d} {mean_len:9.1f} {mean_usable:12.1f} {mean_plddt:11.2f} {collapse*100:9.1f}%")

print()
print("Correlation check (within CONTROL only, where collapse still varies 26/50): does length")
print("predict collapse even absent any steering?")
ctrl = by_cond["CONTROL"]
lens = [float(r["gen_length"]) for r in ctrl]
coll = [float(r["collapse"]) for r in ctrl]
n = len(lens)
mean_l, mean_c = sum(lens) / n, sum(coll) / n
cov = sum((l - mean_l) * (c - mean_c) for l, c in zip(lens, coll)) / n
sd_l = (sum((l - mean_l) ** 2 for l in lens) / n) ** 0.5
sd_c = (sum((c - mean_c) ** 2 for c in coll) / n) ** 0.5
r = cov / (sd_l * sd_c) if sd_l > 0 and sd_c > 0 else float("nan")
print(f"CONTROL point-biserial r(length, collapse), n={n}: {r:+.3f}")

print()
print("Across all conditions pooled:")
all_rows = [r for rows in by_cond.values() for r in rows]
lens = [float(r["gen_length"]) for r in all_rows]
coll = [float(r["collapse"]) for r in all_rows]
n = len(lens)
mean_l, mean_c = sum(lens) / n, sum(coll) / n
cov = sum((l - mean_l) * (c - mean_c) for l, c in zip(lens, coll)) / n
sd_l = (sum((l - mean_l) ** 2 for l in lens) / n) ** 0.5
sd_c = (sum((c - mean_c) ** 2 for c in coll) / n) ** 0.5
r_all = cov / (sd_l * sd_c) if sd_l > 0 and sd_c > 0 else float("nan")
print(f"ALL CONDITIONS pooled point-biserial r(length, collapse), n={n}: {r_all:+.3f}")
print("(Pooled r conflates 'steering causes both more length and more collapse' with any direct")
print(" length effect -- the CONTROL-only number above is the cleaner read on whether length")
print(" alone predicts folding failure.)")
