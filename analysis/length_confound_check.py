"""
Why did ProtGPT2 come back "INERT" in NB49 when SS1g/SS1l saw catastrophic collapse at the
identical absolute push (1167.996)?

Answer: it did NOT come back inert. Steering cut ProtGPT2's output length by ~60%, and pLDDT is
strongly length-dependent -- short peptides score HIGH. So the collapse metric
`collapse = int(0.0 < plddt < 60.0)` scored the damaged output as a success.

This is the mirror image of the SS1l ORTHOGONAL_2x artifact. There, steering made generations
unusually LONG, ESMFold OOM'd, plddt came back 0.0, and the metric scored a total fold failure as
NOT-collapsed. CLAUDE.md sec 5 warns about that direction only. This is the SHORT direction, and
nothing in the project was watching for it.

Zero-GPU. Requires scipy (present locally); no pandas.
Run:  PYTHONIOENCODING=utf-8 python analysis/length_confound_check.py
"""

import csv
import collections
import math
import os

from scipy.stats import fisher_exact, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = ["protgpt2", "zymctrl", "rita", "mistralprot", "progen2"]


def load(model):
    with open(os.path.join(HERE, f"{model}_layer_matched_sequences.csv"),
              newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["L"] = int(r["usable_length"])
        r["p"] = float(r["plddt"])
        r["c"] = int(r["collapse"])
    return rows


def spearman(xs, ys):
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            a = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[o[k]] = a
            i = j + 1
        return r
    rx, ry = rk(xs), rk(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def hr(c="=", n=94):
    print(c * n)


hr()
print("1. DOES STEERING CHANGE GENERATION LENGTH? (mean usable residues, relative to CONTROL)")
hr()
print(f"{'model':<13}{'ctrl len':>10}{'min ratio':>11}{'max ratio':>11}{'Spearman(len,pLDDT)':>22}")
print("-" * 94)
shift = {}
for m in MODELS:
    rows = load(m)
    by = collections.defaultdict(list)
    for r in rows:
        by[r["condition"]].append(r)
    cl = sum(r["L"] for r in by["CONTROL"]) / len(by["CONTROL"])
    ratios = [(sum(r["L"] for r in g) / len(g)) / cl
              for c, g in by.items() if c != "CONTROL"]
    sp = spearman([r["L"] for r in rows], [r["p"] for r in rows])
    shift[m] = (cl, min(ratios), max(ratios), sp)
    print(f"{m:<13}{cl:10.1f}{min(ratios):11.2f}{max(ratios):11.2f}{sp:+22.3f}")
print()
print("ProtGPT2 loses ~60% of its length under steering; Mistral-Prot roughly doubles.")
print("ZymCTRL and RITA do not move at all -- their nulls are clean and unaffected by this.")

hr()
print("2. IS pLDDT ITSELF LENGTH-DRIVEN? (ProtGPT2, all 250 sequences, ignoring condition)")
hr()
rows = load("protgpt2")
print(f"{'length bin':<14}{'n':>5}{'mean pLDDT':>12}{'collapse':>11}")
print("-" * 94)
for lo, hi in [(0, 25), (25, 40), (40, 60), (60, 100), (100, 10 ** 9)]:
    g = [r for r in rows if lo <= r["L"] < hi]
    if not g:
        continue
    lab = f"{lo}-{hi if hi < 10 ** 9 else '+'}"
    print(f"{lab:<14}{len(g):>5}{sum(r['p'] for r in g) / len(g):12.2f}"
          f"{sum(r['c'] for r in g) / len(g) * 100:10.1f}%")

sh = [r for r in rows if r["L"] < 30]
lg = [r for r in rows if r["L"] >= 30]
ks, kl = sum(r["c"] for r in sh), sum(r["c"] for r in lg)
_, pf = fisher_exact([[ks, len(sh) - ks], [kl, len(lg) - kl]])
print()
print(f"  <30 residues : n={len(sh):3d}  pLDDT {sum(r['p'] for r in sh)/len(sh):.2f}  "
      f"collapse {ks/len(sh):.1%}")
print(f"  >=30 residues: n={len(lg):3d}  pLDDT {sum(r['p'] for r in lg)/len(lg):.2f}  "
      f"collapse {kl/len(lg):.1%}")
print(f"  Fisher p = {pf:.4g}")
print("  Sequence length alone predicts 'collapse', with no reference to any condition.")

hr()
print("3. THE VERDICT REVERSES ONCE LENGTH IS HELD ROUGHLY FIXED (ProtGPT2)")
hr()
ctrl = [r for r in rows if r["condition"] == "CONTROL"]
steer = [r for r in rows if r["condition"] != "CONTROL"]


def compare(a, b, label):
    ka, kb = sum(r["c"] for r in a), sum(r["c"] for r in b)
    _, pf = fisher_exact([[ka, len(a) - ka], [kb, len(b) - kb]])
    _, pm = mannwhitneyu([r["p"] for r in a], [r["p"] for r in b], alternative="two-sided")
    print(f"  {label}")
    print(f"    steered  n={len(a):3d}  collapse {ka/len(a):6.1%}  pLDDT {sum(r['p'] for r in a)/len(a):6.2f}")
    print(f"    control  n={len(b):3d}  collapse {kb/len(b):6.1%}  pLDDT {sum(r['p'] for r in b)/len(b):6.2f}")
    print(f"    Fisher p = {pf:.4g}   MWU p = {pm:.4g}")
    print()


compare(steer, ctrl, "ALL lengths  <-- this is what NB49 reported, and why it said INERT")
compare([r for r in steer if r["L"] >= 30], [r for r in ctrl if r["L"] >= 30],
        "length-matched, usable >= 30 residues")
compare([r for r in steer if r["L"] >= 40], [r for r in ctrl if r["L"] >= 40],
        "length-matched, usable >= 40 residues")

hr()
print("4. WHAT THIS DOES AND DOES NOT ESTABLISH")
hr()
print("ESTABLISHED:")
print("  * Steering had a large effect on ProtGPT2 -- it cut output length by ~60%.")
print("  * pLDDT is length-dependent in this data, and length alone predicts the collapse label.")
print("  * Among sequences >= 40 residues, steered output loses ~8 pLDDT points vs control")
print("    (MWU p = 0.016). The 'INERT' verdict is an artifact of the length shift.")
print("  * ZymCTRL and RITA show NO length shift, so their nulls stand as reported.")
print()
print("NOT ESTABLISHED -- do not overclaim:")
print("  * This is post-hoc conditioning on a variable steering itself causes. Length is a")
print("    MEDIATOR, not a nuisance covariate, so the length-matched contrast is suggestive")
print("    evidence that the aggregate null is untrustworthy -- NOT a clean causal estimate of")
print("    'how much damage steering does'. The proper fix is to control length BY CONSTRUCTION.")
print("  * n is small in the matched subsets (28 control vs 42 steered at >= 40).")
print("  * Collapse rate alone stays non-significant (p = 0.14); only pLDDT reaches p < 0.05.")
print("  * SS1g never saved sequence-level data, so its length distribution CANNOT be checked.")
print("    Whether SS1g's L12_2x had the same shift is unknowable from what was kept.")
print()
print("CONSEQUENCE FOR THE PROJECT:")
print("  Every steering condition in which generation length moved is affected, old and new.")
print("  The fix for any future run: generate to a FIXED residue count, or report length as a")
print("  primary outcome next to pLDDT, or stratify. Reporting a bare collapse rate while length")
print("  is free to move measures length as much as it measures structure.")
