"""
Screening efficiency of the collapse predictor — how much folding compute does it save?

Zero-GPU. Uses out-of-sample predictions already saved in
analysis/subclinical_detection_sequences.csv: 50 CONTROL (unsteered) ProtGPT2 generations scored by
the layer-30 classifier that was trained on a SEPARATE N=250 pool (notebook 42). So these are
genuine held-out predictions, not in-sample fits.

The question this answers is the practitioner's question, not the ML one:
    "I can only afford to fold K of my N candidates. If I use the classifier to pick which K,
     how many of the genuinely-foldable ones do I keep — versus just picking at random?"

Reported as enrichment over random, which is the standard way virtual-screening triage is scored
and needs no wall-clock timing to be meaningful.

Run:  python analysis/screening_efficiency.py
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "subclinical_detection_sequences.csv")


def load_control():
    rows = []
    with open(PATH, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["condition"] != "CONTROL":
                continue
            if not r.get("predicted_prob"):
                continue
            rows.append({
                "prob_collapse": float(r["predicted_prob"]),
                "collapse": int(r["collapse"]),
                "plddt": float(r["plddt"]),
            })
    return rows


def hr(c="=", n=76):
    print(c * n)


def main():
    rows = load_control()
    n = len(rows)
    n_good = sum(1 for r in rows if r["collapse"] == 0)
    n_bad = n - n_good
    base_rate = n_good / n

    hr()
    print("SCREENING EFFICIENCY — using the collapse predictor to choose what to fold")
    hr()
    print(f"Pool: {n} unsteered ProtGPT2 generations (held out from the classifier's training set)")
    print(f"Genuinely foldable (pLDDT >= 60): {n_good}/{n} = {base_rate:.1%}")
    print(f"Collapsed:                        {n_bad}/{n}")
    print()

    # Rank ascending by predicted collapse probability = most-promising first.
    ranked = sorted(rows, key=lambda r: r["prob_collapse"])

    print("If you can only afford to fold the top-K most promising candidates:")
    print()
    print(f"{'Budget':>8s} {'Folded':>7s} {'Good found':>11s} {'Recall':>8s} "
          f"{'Precision':>10s} {'Random':>8s} {'Enrich':>7s}")
    print("-" * 68)

    results = []
    for frac in (0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00):
        k = max(1, int(round(frac * n)))
        top = ranked[:k]
        found = sum(1 for r in top if r["collapse"] == 0)
        recall = found / n_good if n_good else float("nan")
        precision = found / k
        expected_random = base_rate * k          # what random selection would find
        enrichment = (found / expected_random) if expected_random > 0 else float("nan")
        results.append((frac, k, found, recall, precision, enrichment))
        print(f"{frac:7.0%} {k:7d} {found:5d}/{n_good:<5d} {recall:7.1%} "
              f"{precision:9.1%} {expected_random:7.1f} {enrichment:6.2f}x")

    print()
    hr("-")
    print("THE COMPUTE QUESTION, ANSWERED TWO WAYS")
    hr("-")

    # (a) Fold-budget needed to retain most of the good candidates.
    for target in (0.80, 0.90, 1.00):
        need = None
        for frac, k, found, recall, precision, enr in results:
            if recall >= target - 1e-9:
                need = (frac, k, recall)
                break
        if need:
            frac, k, recall = need
            saved = 1.0 - frac
            print(f"  Retain {target:.0%} of foldable candidates -> fold only {frac:.0%} "
                  f"({k}/{n}). Folding compute saved: {saved:.0%}")
        else:
            print(f"  Retain {target:.0%}: not achievable below 100% budget in this pool")

    # (b) Yield improvement at a fixed budget.
    print()
    half = [r for r in results if abs(r[0] - 0.50) < 1e-9][0]
    print(f"  At a fixed 50% fold budget: {half[4]:.1%} of what you fold is actually good, "
          f"vs {base_rate:.1%} picking at random ({half[5]:.2f}x enrichment).")

    print()
    hr("-")
    print("HONEST LIMITS OF THIS NUMBER — state these if it goes in the paper")
    hr("-")
    print(f"  * N={n} (single condition, one model). Wide binomial error at every point;")
    print(f"    treat the curve as indicative, not a precision estimate.")
    print("  * Uses FULL-SEQUENCE features. The early-abort result (locked-results.md §3g) uses")
    print("    truncated features and has a LOWER AUC (0.740 at 30 residues vs 0.780 full), so")
    print("    the savings above are an UPPER bound for a true abort-early-mid-generation system.")
    print("  * 'Good' here means ESMFold pLDDT >= 60 — a confidence proxy, not verified folding.")
    print("  * Says nothing about wall-clock cost ratios (generation vs folding vs classifier);")
    print("    that needs a timing measurement, which no notebook has recorded yet.")
    print()


if __name__ == "__main__":
    main()
