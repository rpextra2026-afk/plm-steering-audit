"""
Pre-writing verification: does every headline number the paper will cite actually match its
source artifact?

This is the last check before drafting. It re-derives each key figure from the CSV that produced it
and compares against what notes/locked-results.md claims, so a transcription slip cannot reach the
paper. Zero-GPU, stdlib only.

Numbers whose source is a notebook output cell (not a CSV) are listed as MANUAL — they were read
from saved cells and cannot be re-derived here; they are flagged so they get eyeballed once.

Run:  PYTHONIOENCODING=utf-8 python analysis/verify_paper_numbers.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NOTES = os.path.join(HERE, "..", "notes", "locked-results.md")

PASS, FAIL, MANUAL = [], [], []


def load(name):
    with open(os.path.join(HERE, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def check(label, got, want, tol, unit=""):
    ok = got is not None and abs(got - want) <= tol
    (PASS if ok else FAIL).append(label)
    flag = "PASS" if ok else "**FAIL**"
    gs = f"{got:.4f}" if got is not None else "n/a"
    print(f"  [{flag}] {label}: computed {gs}{unit}, claimed {want}{unit} (tol {tol})")


def manual(label, value, source):
    MANUAL.append(label)
    print(f"  [MANUAL] {label} = {value}  <- from {source}, not re-derivable from CSV")


def hr(c="=", n=86):
    print(c * n)


def main():
    notes = open(NOTES, encoding="utf-8").read()

    hr()
    print("VERIFYING PAPER NUMBERS AGAINST SOURCE ARTIFACTS")
    hr()

    # ---------------------------------------------------------------- compute cost (§3k)
    print("\n§3k — measured compute cost (source: compute_cost_summary.csv)")
    cs = load("compute_cost_summary.csv")[0]
    t_gen, t_fold = float(cs["t_gen_full_s"]), float(cs["t_fold_s"])
    check("ESMFold seconds/sequence", t_fold, 2.095, 0.001, "s")
    check("Generation seconds/sequence", t_gen, 0.592, 0.001, "s")
    check("fold/generation ratio", t_fold / t_gen, 3.54, 0.01, "x")
    check("fold share of pipeline", t_fold / (t_fold + t_gen), 0.78, 0.005)
    check("fold vs classifier ratio", t_fold / float(cs["t_clf_30_s"]), 77.0, 1.0, "x")
    check("natural collapse rate", float(cs["collapse_rate"]), 0.50, 0.005)
    check("AUC full (pooled OOF)", float(cs["auc_full"]), 0.724, 0.001)
    check("AUC 30 residues", float(cs["auc_30"]), 0.715, 0.001)

    # ---------------------------------------------------------------- triage (§3k)
    print("\n§3k — unbatched triage policies (source: compute_triage_policies.csv)")
    tp = load("compute_triage_policies.csv")
    row = [r for r in tp if int(float(r["abort_at"])) == 30 and abs(float(r["budget"]) - 0.7) < 1e-6]
    if row:
        check("abort@30 b=70% recall", float(row[0]["recall"]), 0.808, 0.002)
        check("abort@30 b=70% saved", float(row[0]["frac_saved"]), 0.248, 0.002)
    else:
        FAIL.append("abort@30 b=70% row missing")
        print("  [**FAIL**] abort@30 b=70% row not found in compute_triage_policies.csv")

    # ---------------------------------------------------------------- batched (§3k-BATCHED)
    print("\n§3k-BATCHED — batched throughput (source: batched_*.csv)")
    gen = {int(r["batch"]): r for r in load("batched_generation_throughput.csv")}
    fold = load("batched_folding_throughput.csv")
    sortedf = {int(r["batch"]): r for r in fold if r["ordering"] == "length-sorted"}
    naivef = {int(r["batch"]): r for r in fold if r["ordering"] == "naive order"}

    check("generation speedup @batch16", float(gen[16]["speedup"]), 5.80, 0.02, "x")
    check("folding speedup @batch4 sorted", float(sortedf[4]["speedup"]), 1.33, 0.01, "x")
    check("naive batch4 is SLOWER than batch1", float(naivef[4]["speedup"]), 0.755, 0.01, "x")
    check("peak VRAM folding @batch4", float(sortedf[4]["peak_gb"]), 15.06, 0.02, "GB")

    tg_b, tf_b = float(gen[16]["s_per_seq"]), float(sortedf[4]["s_per_seq"])
    check("fold share, both batched", tf_b / (tg_b + tf_b), 0.95, 0.005)

    # independently re-derive the headline saving from raw timings
    T_CLF = float(cs["t_clf_30_s"])
    frac_gen = 30.0 / 47.3
    tg_k = tg_b * frac_gen
    base = tg_b + tf_b
    cost70 = tg_k + T_CLF + 0.70 * ((tg_b - tg_k) + tf_b)
    check("HEADLINE saving @b=70% (re-derived)", 1 - cost70 / base, 0.273, 0.005)

    # worst-case bound quoted in the notes
    tg_unb = float(gen[1]["s_per_seq"])
    cost70_worst = tg_k + T_CLF + 0.70 * (tg_unb * (1 - frac_gen) + tf_b)
    check("worst-case bound @b=70%", 1 - cost70_worst / base, 0.211, 0.006)

    # ---------------------------------------------------------------- validity (§3k-BATCHED)
    print("\n§3k-BATCHED — batching validity (source: batched_validity_check.csv)")
    v = load("batched_validity_check.csv")
    flips = sum(int(r["label_flips"]) for r in v)
    check("total collapse-label flips", float(flips), 0.0, 0.0)
    check("max abs pLDDT diff", max(float(r["max_abs_diff"]) for r in v), 0.051, 0.001)

    # ---------------------------------------------------------------- early abort (§3j)
    print("\n§3j — early abort across models (source: early_abort_multimodel_*.csv)")
    comb = {r["model"]: r for r in load("early_abort_multimodel_combined.csv")}
    zy = comb.get("ZymCTRL")
    ri = comb.get("RITA-small")
    if zy:
        check("ZymCTRL AUC full", float(zy["auc_full"]), 0.658, 0.002)
        check("ZymCTRL AUC 30res", float(zy["auc_30"]), 0.646, 0.002)
        check("ZymCTRL retention@30", float(zy["retention_30"]), 0.981, 0.002)
    if ri:
        check("RITA AUC full (must be ~chance)", float(ri["auc_full"]), 0.506, 0.002)
        # the critical assertion: RITA must NOT be counted as a working model
        ok = abs(float(ri["auc_full"]) - 0.5) < 0.05
        (PASS if ok else FAIL).append("RITA correctly at chance")
        print(f"  [{'PASS' if ok else '**FAIL**'}] RITA is at chance (|AUC-0.5| < 0.05) "
              f"-> must be reported as a NULL, not a third success")

    # ---------------------------------------------------------------- screening (preliminary)
    print("\nPreliminary screening (source: subclinical_detection_sequences.csv)")
    rows = [r for r in load("subclinical_detection_sequences.csv") if r["condition"] == "CONTROL"]
    n_good = sum(1 for r in rows if r["collapse"] == "0")
    check("CONTROL pool size", float(len(rows)), 50.0, 0.0)
    check("CONTROL foldable count", float(n_good), 19.0, 0.0)

    # ---------------------------------------------------------------- notes consistency
    print("\nCross-checking claims present in locked-results.md text")
    for needle, label in [
        ("3.54", "fold/gen ratio 3.54x appears in notes"),
        ("77", "77x classifier ratio appears in notes"),
        ("27.3%", "headline batched saving 27.3% appears in notes"),
        ("0/64", "zero label flips appears in notes"),
        ("95.8%", "top-decile precision 95.8% appears in notes"),
        ("1.24", "dose50 = 1.24x appears in notes"),
    ]:
        ok = needle in notes
        (PASS if ok else FAIL).append(label)
        print(f"  [{'PASS' if ok else '**FAIL**'}] {label}")

    # ---------------------------------------------------------------- manual-only figures
    print("\nFigures read from notebook output cells (NOT re-derivable here — eyeball once)")
    manual("§3a raw layer-30 AUC = 0.772 / entropy 0.586 / OSAE 0.528", "0.772", "NB08 cell output")
    manual("§3c top-decile precision = 95.8%", "95.8%", "NB08 appended cells")
    manual("§3g ProtGPT2 early-abort 0.682/0.720/0.740/0.780", "see §3g", "NB42 cell output")
    manual("§1n dose50 = 1.240x, CI [1.138, 1.321]", "1.240x", "NB37 cell output")
    manual("§1l random/shuffled/orthogonal all 100% at 2x", "100%", "NB35/NB39 cell outputs")

    hr()
    print("SUMMARY")
    hr()
    print(f"  PASS   : {len(PASS)}")
    print(f"  FAIL   : {len(FAIL)}")
    print(f"  MANUAL : {len(MANUAL)} (read from notebook cells; verify by eye once)")
    if FAIL:
        print("\n  FAILURES:")
        for f in FAIL:
            print(f"    - {f}")
        print("\n  ==> DO NOT DRAFT until these are reconciled.")
    else:
        print("\n  ==> All CSV-derivable numbers match the notes. Safe to draft from")
        print("      notes/locked-results.md.")


if __name__ == "__main__":
    main()
