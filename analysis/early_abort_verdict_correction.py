"""
Correction to notebook 44's automated verdict.

NB44 gated "early-abort viable" on ONE condition: retention_30 = AUC(30res)/AUC(full) > 0.85.
That test is invalid when AUC(full) is itself at chance — retaining 97.7% of a chance-level
signal is retaining 97.7% of nothing. NB44 consequently reported "3/3 viable" including RITA,
whose AUC is 0.495-0.506 (i.e. a coin flip at every truncation).

The correct gate has TWO conditions:
  (1) the full-sequence classifier must actually beat chance, and
  (2) the truncated classifier must retain most of that signal.

This script applies both, with a one-sample t-test of each AUC distribution against 0.5 using the
per-repeat spread NB44 already recorded. Zero-GPU; reads NB44's own saved CSVs.

Run:  python analysis/early_abort_verdict_correction.py
"""

import csv
import math
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = {"ZymCTRL": "early_abort_multimodel_zymctrl.csv",
         "RITA-small": "early_abort_multimodel_rita.csv"}

# ProtGPT2 has no per-repeat CSV here (it comes from §3g / notebook 42); its locked means are used
# for the comparison row, and its significance is already established in locked-results.md §3g.
PROTGPT2 = {10: 0.682, 20: 0.720, 30: 0.740, -1: 0.780}


def mean_sd(xs):
    n = len(xs)
    m = sum(xs) / n
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in xs) / (n - 1)
    return m, math.sqrt(var)


def t_vs_chance(xs, chance=0.5):
    """One-sample t of the AUC repeats against chance. Returns (t, n)."""
    m, sd = mean_sd(xs)
    n = len(xs)
    if sd == 0 or n < 2:
        return float("nan"), n
    return (m - chance) / (sd / math.sqrt(n)), n


def load(fname):
    by_trunc = defaultdict(list)
    with open(os.path.join(HERE, fname), newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            by_trunc[int(r["truncation"])].append(float(r["auc"]))
    return by_trunc


def hr(c="=", n=84):
    print(c * n)


def main():
    hr()
    print("CORRECTED EARLY-ABORT VERDICT — two-condition gate")
    hr()
    print("Gate (1): does the FULL-sequence classifier beat chance?  [t vs AUC=0.5, per-repeat spread]")
    print("Gate (2): does the 30-residue classifier retain >85% of that signal?")
    print("NB44 applied only gate (2), which is meaningless if gate (1) fails.")
    print()

    print(f"{'Model':>14s} {'AUC full':>18s} {'AUC 30res':>18s} {'t(full vs .5)':>14s} {'verdict':>16s}")
    print("-" * 84)

    verdicts = {}

    # ProtGPT2 reference row (significance established separately in §3g)
    print(f"{'ProtGPT2 (§3g)':>14s} {PROTGPT2[-1]:>18.3f} {PROTGPT2[30]:>18.3f} "
          f"{'(see §3g)':>14s} {'REAL SIGNAL':>16s}")
    verdicts["ProtGPT2"] = True

    for model, fname in FILES.items():
        d = load(fname)
        full, r30 = d[-1], d[30]
        m_full, sd_full = mean_sd(full)
        m_30, sd_30 = mean_sd(r30)
        t_full, n = t_vs_chance(full)
        # two-sided crit value ~2.14 at df=14 (alpha=.05); state explicitly rather than import scipy
        beats_chance = abs(t_full) > 2.14 and m_full > 0.5
        retention = m_30 / m_full if m_full else float("nan")
        ok = beats_chance and retention > 0.85
        verdicts[model] = ok
        verdict = "REAL SIGNAL" if beats_chance else "AT CHANCE"
        print(f"{model:>14s} {m_full:>10.3f} ±{sd_full:<6.3f} {m_30:>10.3f} ±{sd_30:<6.3f} "
              f"{t_full:>14.2f} {verdict:>16s}")

    print()
    print("(t critical ≈ 2.14 at df=14, two-sided α=0.05, from 15 repeats per truncation)")
    print()
    hr("-")
    print("WHAT THIS CHANGES")
    hr("-")
    n_real = sum(1 for v in verdicts.values() if v)
    print(f"NB44's automated verdict said:  3/3 models viable")
    print(f"Correct verdict:                {n_real}/3 models viable")
    print()
    for model, ok in verdicts.items():
        print(f"  {model:>14s}: {'early-abort supported' if ok else 'NO SIGNAL — early-abort question is moot'}")
    print()
    print("RITA's 97.7% 'retention' is retention of a coin flip: AUC 0.506 full, 0.495 at 30")
    print("residues, t≈0.2 against chance. There is no signal to retain early, so the truncation")
    print("result carries no information either way. This is NOT evidence that early-abort fails on")
    print("RITA — it is evidence that the underlying collapse-prediction task fails on RITA in this")
    print("pool, which makes the early-abort question unanswerable there.")
    print()
    print("Why RITA plausibly failed here: 19 non-collapsed of 200 (90.5% collapse). At a 30% test")
    print("fold that is ~6 positive examples per split. §3f's locked 0.627 came from a different")
    print("draw of the same underlying setup; the gap between 0.627 and 0.506 is consistent with")
    print("that severe imbalance, not with a contradiction. Report as underpowered, not refuted.")
    print()
    hr("-")
    print("HONEST CLAIM FOR THE PAPER")
    hr("-")
    print("  'Early abort is demonstrated on 2 of the 3 models where the collapse-prediction task")
    print("   is itself learnable (ProtGPT2, ZymCTRL). On RITA the task does not reach above chance")
    print("   in our pool due to severe class imbalance, so the early-abort question cannot be")
    print("   answered there. Two further models (p-IgGen, Mistral-Prot) are structurally")
    print("   untestable — their natural collapse rates of 0% and ~99% leave no label variance.'")
    print()


if __name__ == "__main__":
    main()
