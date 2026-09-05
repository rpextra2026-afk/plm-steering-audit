"""
Local label sanity-check (pre-writing-lock-checklist item 12; context-and-decisions.md §6).

Reads already-saved generations with their ESMFold pLDDT + collapse label and prints the
most-collapsed vs. most-confident CONTROL (unsteered) sequences side by side, with cheap
degeneracy diagnostics (longest homopolymer run, distinct-3-mer fraction, Shannon entropy
already in the file), to confirm the pLDDT<60 collapse LABEL matches visible sequence reality
before any of it goes in the paper. This is the "read the actual sequences" gut-check a
reviewer may ask for.

No GPU, no model, no numpy — pure stdlib over analysis/*_sequences.csv.
Run:  python analysis/eyeball_sequences.py
"""

import csv
from collections import Counter

SRC = "analysis/orthogonalized_foldability_and_refold_sequences.csv"
COLLAPSE_THRESHOLD = 60.0  # pLDDT < 60 = "collapsed", the project's primary threshold


def longest_homopolymer(s):
    best = cur = 0
    prev = None
    for c in s:
        cur = cur + 1 if c == prev else 1
        prev = c
        best = max(best, cur)
    return best


def distinct_kmer_frac(s, k=3):
    if len(s) < k:
        return float("nan")
    kmers = [s[i:i + k] for i in range(len(s) - k + 1)]
    return len(set(kmers)) / len(kmers)


rows = []
with open(SRC, newline="") as f:
    for r in csv.DictReader(f):
        if r.get("condition") != "CONTROL":
            continue  # CONTROL = unsteered natural ProtGPT2 generation
        raw = r.get("gen_only") or r.get("sequence") or ""
        seq = "".join(raw.split())  # strip the CSV's line-wrap whitespace
        try:
            plddt = float(r["plddt"])
        except (KeyError, ValueError):
            continue
        try:
            entropy = float(r.get("entropy", ""))
        except ValueError:
            entropy = float("nan")
        rows.append({
            "plddt": plddt,
            "collapse": r.get("collapse"),
            "entropy": entropy,
            "seq": seq,
            "hp": longest_homopolymer(seq),
            "d3": distinct_kmer_frac(seq),
            "len": len(seq),
        })

rows.sort(key=lambda x: x["plddt"])
n_collapsed = sum(1 for r in rows if r["plddt"] < COLLAPSE_THRESHOLD)
print(f"SRC: {SRC}")
print(f"CONTROL sequences with a pLDDT label: {len(rows)}")
print(f"  collapsed (pLDDT<{COLLAPSE_THRESHOLD:.0f}): {n_collapsed}"
      f"   healthy (>={COLLAPSE_THRESHOLD:.0f}): {len(rows) - n_collapsed}")


def show(title, subset):
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)
    for r in subset:
        print(f"pLDDT={r['plddt']:5.1f}  collapse={r['collapse']}  entropy={r['entropy']:.2f}"
              f"  max_homopolymer={r['hp']}  distinct_3mer={r['d3']:.2f}  len={r['len']}")
        print("   " + r["seq"])


show("5 LOWEST pLDDT  (label = collapsed)", rows[:5])
show("5 HIGHEST pLDDT (label = healthy)", rows[-5:])

# Group-level contrast: does the label track the cheap degeneracy proxies at all?
def mean(xs):
    xs = [x for x in xs if x == x]  # drop NaN
    return sum(xs) / len(xs) if xs else float("nan")

lo = [r for r in rows if r["plddt"] < COLLAPSE_THRESHOLD]
hi = [r for r in rows if r["plddt"] >= COLLAPSE_THRESHOLD]
print("\n" + "-" * 72)
print("group means (collapsed vs healthy):")
print(f"  entropy:        {mean([r['entropy'] for r in lo]):.3f}  vs  {mean([r['entropy'] for r in hi]):.3f}")
print(f"  max_homopolymer:{mean([r['hp'] for r in lo]):.3f}  vs  {mean([r['hp'] for r in hi]):.3f}")
print(f"  distinct_3mer:  {mean([r['d3'] for r in lo]):.3f}  vs  {mean([r['d3'] for r in hi]):.3f}")
