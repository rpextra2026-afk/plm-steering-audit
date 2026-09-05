"""
Are the predictor's disagreements with ESMFold RANDOM or SYSTEMATIC?

Context. The layer-30 classifier is trained on pLDDT labels, so it is downstream of ESMFold, not a
competitor to it — it cannot beat the teacher it is graded against. Before investing in real
experimental labels to break that circularity, this asks a cheaper question:

    When the classifier disagrees with ESMFold, is the disagreement structured?

  * RANDOM disagreement  -> the classifier is a noisy pLDDT imitator. Nothing to chase. Idea dies
                            cheap, which is a real result.
  * SYSTEMATIC disagreement -> the PLM's internal state tracks something pLDDT does not, and the
                            "independent confidence score" direction is worth real data.

This CANNOT show the classifier is *better* than pLDDT — that needs experimental ground truth
(stability/expression), by construction. It only screens for whether anything is there.

Data: analysis/subclinical_detection_sequences.csv — 50 CONTROL (unsteered) + 50 REAL_1x ProtGPT2
generations, each with a HELD-OUT predicted probability (classifier trained on a separate 250-pool,
notebook 42) and its real ESMFold pLDDT.

Zero-GPU, stdlib only.
Run:  PYTHONIOENCODING=utf-8 python analysis/plddt_disagreement_analysis.py
"""

import csv
import math
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "subclinical_detection_sequences.csv")

HYDROPHOBIC = set("AVILMFWY")
CHARGED = set("DEKR")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------- sequence features
def shannon(seq):
    if not seq:
        return 0.0
    c = Counter(seq); n = len(seq)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def distinct_kmer(seq, k=3):
    if len(seq) < k:
        return 1.0
    grams = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    return len(set(grams)) / len(grams)


def longest_homopolymer(seq):
    if not seq:
        return 0
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def frac(seq, group):
    return sum(1 for a in seq if a in group) / len(seq) if seq else 0.0


def features(seq):
    s = "".join(a for a in seq if a in VALID_AA)
    return {
        "length": len(s),
        "entropy": shannon(s),
        "distinct3": distinct_kmer(s, 3),
        "max_homopolymer": longest_homopolymer(s),
        "hydrophobic_frac": frac(s, HYDROPHOBIC),
        "charged_frac": frac(s, CHARGED),
    }


# ---------------------------------------------------------------- stats
def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def sd(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = mean(xs), mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return num / den if den else float("nan")


def t_from_r(r, n):
    if abs(r) >= 1 or n < 3:
        return float("nan")
    return r * math.sqrt((n - 2) / (1 - r * r))


def welch(a, b):
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    va, vb = sd(a) ** 2 / len(a), sd(b) ** 2 / len(b)
    return (mean(a) - mean(b)) / math.sqrt(va + vb) if (va + vb) > 0 else float("nan")


def hr(c="=", n=88):
    print(c * n)


def analyse(rows, label):
    hr()
    print(f"DISAGREEMENT ANALYSIS — {label}  (n={len(rows)})")
    hr()

    # pLDDT rescaled to a 0-1 "ESMFold says it folds" score, so both scores point the same way.
    # predicted_prob is P(collapse), so 1-p is "classifier says it folds".
    for r in rows:
        r["clf_fold"] = 1.0 - r["prob"]
        r["esm_fold"] = min(1.0, max(0.0, r["plddt"] / 100.0))

    cf = [r["clf_fold"] for r in rows]
    ef = [r["esm_fold"] for r in rows]
    r_agree = pearson(cf, ef)
    print(f"Agreement between the two scores: r = {r_agree:+.3f} "
          f"(t={t_from_r(r_agree, len(rows)):+.2f}, n={len(rows)})")
    print("  (a classifier that perfectly imitated ESMFold would approach r=+1)")

    # Residual: how much MORE optimistic the classifier is than ESMFold, after removing the
    # shared linear trend. Positive = classifier says fold, ESMFold disagrees.
    mx, my = mean(ef), mean(cf)
    den = sum((x - mx) ** 2 for x in ef)
    slope = sum((x - mx) * (y - my) for x, y in zip(ef, cf)) / den if den else 0.0
    icept = my - slope * mx
    for r in rows:
        r["resid"] = r["clf_fold"] - (icept + slope * r["esm_fold"])

    resids = [r["resid"] for r in rows]
    print(f"Residual spread: sd = {sd(resids):.3f}  (0 would mean perfect imitation)")

    print()
    print("Is the disagreement STRUCTURED? — correlation of residual with sequence features")
    print(f"{'feature':>18s} {'r':>8s} {'t':>8s} {'verdict':>14s}")
    print("-" * 54)
    feat_names = ["length", "entropy", "distinct3", "max_homopolymer",
                  "hydrophobic_frac", "charged_frac"]
    hits = []
    for f in feat_names:
        xs = [r["feat"][f] for r in rows]
        r_f = pearson(xs, resids)
        t_f = t_from_r(r_f, len(rows))
        sig = abs(t_f) > 2.0
        if sig:
            hits.append((f, r_f, t_f))
        print(f"{f:>18s} {r_f:+8.3f} {t_f:+8.2f} {'SYSTEMATIC' if sig else 'ns':>14s}")

    print()
    print("The extreme disagreements (what the two scores fight about most):")
    ranked = sorted(rows, key=lambda r: r["resid"])
    print("\n  Classifier says COLLAPSE, ESMFold says folds (most negative residual):")
    for r in ranked[:3]:
        print(f"    pLDDT {r['plddt']:5.1f} | P(collapse) {r['prob']:.2f} | "
              f"len {r['feat']['length']:3d} ent {r['feat']['entropy']:.2f} "
              f"d3 {r['feat']['distinct3']:.2f} | {r['seq'][:44]}")
    print("\n  Classifier says FOLDS, ESMFold says collapse (most positive residual):")
    for r in ranked[-3:]:
        print(f"    pLDDT {r['plddt']:5.1f} | P(collapse) {r['prob']:.2f} | "
              f"len {r['feat']['length']:3d} ent {r['feat']['entropy']:.2f} "
              f"d3 {r['feat']['distinct3']:.2f} | {r['seq'][:44]}")

    # Contrast the two disagreement camps directly.
    n_tail = max(5, len(rows) // 5)
    low, high = ranked[:n_tail], ranked[-n_tail:]
    print()
    print(f"Comparing the two disagreement tails (n={n_tail} each), Welch t:")
    print(f"{'feature':>18s} {'clf-pessimistic':>17s} {'clf-optimistic':>16s} {'t':>8s}")
    print("-" * 62)
    tail_hits = []
    for f in feat_names:
        a = [r["feat"][f] for r in low]
        b = [r["feat"][f] for r in high]
        t = welch(b, a)
        if abs(t) > 2.0:
            tail_hits.append((f, t))
        print(f"{f:>18s} {mean(a):17.3f} {mean(b):16.3f} {t:+8.2f}")

    return hits, tail_hits, r_agree


def main():
    raw = []
    with open(PATH, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if not r.get("predicted_prob"):
                continue
            seq = r["sequence"]
            raw.append({"cond": r["condition"], "seq": seq,
                        "plddt": float(r["plddt"]), "prob": float(r["predicted_prob"]),
                        "collapse": int(r["collapse"]), "feat": features(seq)})

    control = [r for r in raw if r["cond"] == "CONTROL"]

    hits_c, tails_c, r_c = analyse(control, "CONTROL only (unsteered — the clean set)")
    print()
    hits_a, tails_a, r_a = analyse(raw, "CONTROL + REAL_1x pooled (more power, mixed distribution)")

    print()
    hr()
    print("VERDICT — is the 'independent confidence score' idea worth pursuing?")
    hr()
    all_hits = set(f for f, _, _ in hits_c) | set(f for f, _, _ in hits_a)
    all_tails = set(f for f, _ in tails_c) | set(f for f, _ in tails_a)
    flagged = all_hits | all_tails

    print(f"Score agreement (CONTROL): r = {r_c:+.3f}")
    print(f"Features where disagreement is systematic: "
          f"{', '.join(sorted(flagged)) if flagged else 'NONE'}")
    print()
    if not flagged:
        print("  ==> DISAGREEMENT LOOKS RANDOM. On this data the classifier reads as a noisy")
        print("      imitator of pLDDT — its errors do not track any sequence property tested.")
        print("      That is a genuine (cheap) negative for the 'independent confidence score'")
        print("      idea AS CURRENTLY TRAINED, and it is exactly what circularity predicts:")
        print("      a model trained on pLDDT learns pLDDT, including its blind spots.")
        print()
        print("      IMPORTANT — this does NOT kill the research direction. It says the current")
        print("      pLDDT-trained classifier cannot be repurposed as an independent score. A")
        print("      model trained on EXPERIMENTAL labels is a different object and is untested.")
    else:
        print("  ==> DISAGREEMENT IS STRUCTURED. The classifier's errors track sequence properties,")
        print("      meaning the PLM's internal state carries information pLDDT does not use.")
        print("      This is the encouraging outcome: worth pursuing with real experimental labels.")
        print("      Next step: get stability (e.g. Tsuboyama megascale ΔG) or expression data and")
        print("      test whether PLM activations beat pLDDT at predicting the REAL outcome.")
    print()
    print("Hard limits of this screen, regardless of outcome:")
    print("  * n=50 (CONTROL). Underpowered for weak effects; only strong structure would show.")
    print("  * Cannot demonstrate 'better than pLDDT' — the labels ARE pLDDT. Needs real data.")
    print("  * Absence of correlation with these 6 features is not absence of ALL structure.")


if __name__ == "__main__":
    main()
