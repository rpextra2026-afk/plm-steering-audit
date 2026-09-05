"""
Correctness audit of the compute-saving claim (locked-results.md §3k).

The §3k cost model uses a SINGLE MEAN fold time for every candidate. That is only valid if fold
cost is independent of whether a sequence collapses. If collapsed sequences (the ones triage skips)
are systematically CHEAPER to fold than healthy ones, then skipping them saves less than the mean
implies and §3k OVERSTATES the savings. If they are more expensive, §3k understates.

Fold time was measured to correlate with length at r=+0.892, so this is not a hypothetical concern —
it depends entirely on whether length differs by collapse status.

Also audits: the linear-generation assumption, and how much of the claim survives if ESMFold were
batched instead of run one-at-a-time.

Zero-GPU. Reads analysis/compute_timing_per_sequence.csv (NB45's raw per-sequence timings).

Run:  PYTHONIOENCODING=utf-8 python analysis/compute_claim_audit.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(HERE, "compute_timing_per_sequence.csv")


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def sd(xs):
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def welch_t(a, b):
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    va, vb = sd(a) ** 2 / len(a), sd(b) ** 2 / len(b)
    if va + vb == 0:
        return float("nan"), float("nan")
    t = (mean(a) - mean(b)) / math.sqrt(va + vb)
    df = (va + vb) ** 2 / ((va ** 2 / (len(a) - 1)) + (vb ** 2 / (len(b) - 1)))
    return t, df


def pearson(xs, ys):
    n = len(xs)
    mx, my = mean(xs), mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    return num / den if den else float("nan")


def hr(c="=", n=82):
    print(c * n)


def main():
    rows = []
    with open(PATH, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["fold_ok"].strip().lower() not in ("true", "1"):
                continue
            rows.append({
                "gen_t": float(r["gen_time_s"]), "fold_t": float(r["fold_time_s"]),
                "fold_len": int(float(r["fold_len"])), "n_tok": int(float(r["n_tok_out"])),
                "collapse": int(r["collapse"]), "gen_len": len(r["gen_only"]),
            })

    collapsed = [r for r in rows if r["collapse"] == 1]
    healthy = [r for r in rows if r["collapse"] == 0]

    hr()
    print("AUDIT 1 — THE ONE THAT MATTERS: is fold cost independent of collapse status?")
    hr()
    print(f"Pool: {len(rows)} folded  |  collapsed {len(collapsed)}  |  healthy {len(healthy)}")
    print()
    ct = [r["fold_t"] for r in collapsed]
    ht = [r["fold_t"] for r in healthy]
    cl = [r["fold_len"] for r in collapsed]
    hl = [r["fold_len"] for r in healthy]

    print(f"{'':>12s} {'fold time (s)':>20s} {'folded length':>20s}")
    print(f"{'collapsed':>12s} {mean(ct):>12.3f} ±{sd(ct):<6.3f} {mean(cl):>12.1f} ±{sd(cl):<6.1f}")
    print(f"{'healthy':>12s} {mean(ht):>12.3f} ±{sd(ht):<6.3f} {mean(hl):>12.1f} ±{sd(hl):<6.1f}")

    t_time, df_time = welch_t(ct, ht)
    t_len, df_len = welch_t(cl, hl)
    print()
    print(f"Welch t (fold time, collapsed vs healthy): t={t_time:+.2f}, df~{df_time:.0f}")
    print(f"Welch t (length,    collapsed vs healthy): t={t_len:+.2f}, df~{df_len:.0f}")
    print(f"corr(length, fold time) overall: r={pearson([r['fold_len'] for r in rows], [r['fold_t'] for r in rows]):+.3f}")

    ratio = mean(ct) / mean(ht) if mean(ht) else float("nan")
    print()
    sig = abs(t_time) > 2.0
    if not sig:
        print(f"  ==> NO significant difference in fold cost by collapse status (|t|={abs(t_time):.2f} < 2).")
        print(f"      The single-mean assumption in §3k is SOUND. Savings estimate is unbiased on")
        print(f"      this axis.")
    else:
        direction = "CHEAPER" if mean(ct) < mean(ht) else "MORE EXPENSIVE"
        print(f"  ==> Collapsed sequences are significantly {direction} to fold "
              f"({mean(ct):.3f}s vs {mean(ht):.3f}s, ratio {ratio:.2f}x).")
        if mean(ct) < mean(ht):
            print(f"      This means §3k OVERSTATES savings: triage skips the cheap folds, so the")
            print(f"      real saving is LOWER than the mean-based model implies.")
            print(f"      Corrected upper bound on fold-compute saved when skipping fraction f:")
            print(f"        f x {ratio:.2f} of mean fold cost, not f x 1.00.")
        else:
            print(f"      §3k UNDERSTATES savings: triage skips the expensive folds, so the real")
            print(f"      saving is HIGHER than reported. §3k is conservative — safe to cite.")

    hr()
    print("AUDIT 2 — is the linear-generation assumption OK?")
    hr()
    gl = [r["gen_len"] for r in rows]
    gt = [r["gen_t"] for r in rows]
    r_gen = pearson(gl, gt)
    print(f"corr(generated residues, generation time) = {r_gen:+.3f}")
    print(f"generation time: mean {mean(gt):.4f}s ±{sd(gt):.4f}")
    print(f"generated residues: mean {mean(gl):.1f} ±{sd(gl):.1f}")
    # implied fixed overhead via simple least-squares intercept
    mx, my = mean(gl), mean(gt)
    denom = sum((x - mx) ** 2 for x in gl)
    slope = sum((x - mx) * (y - my) for x, y in zip(gl, gt)) / denom if denom else float("nan")
    intercept = my - slope * mx
    print(f"least-squares fit: time ~ {intercept:.4f}s + {slope:.5f}s/residue")
    print(f"  fixed overhead is {intercept/my:.0%} of mean generation time")
    if intercept > 0.25 * my:
        print(f"  ==> Non-trivial fixed overhead. §3k's T_gen(K) = per_token x K UNDERSTATES the cost")
        print(f"      of a short generation, which makes early-abort look slightly BETTER than it is.")
        print(f"      Magnitude is small (generation is only 22% of pipeline), but state the")
        print(f"      assumption rather than leaving it implicit.")
    else:
        print(f"  ==> Fixed overhead is small; linear scaling is a fair approximation.")

    hr()
    print("AUDIT 3 — how much of the claim survives if ESMFold were BATCHED?")
    hr()
    t_fold, t_gen = mean([r["fold_t"] for r in rows]), mean(gt)
    print(f"Measured one-at-a-time: fold {t_fold:.3f}s ({t_fold/(t_fold+t_gen):.0%} of pipeline)")
    print()
    print(f"{'batch speedup':>14s} {'fold share':>12s} {'max saving @70% budget':>24s}")
    for sp in (1, 2, 4, 8):
        tf = t_fold / sp
        share = tf / (tf + t_gen)
        # best case: skip 30% of folds entirely, generation unchanged
        saving = 0.30 * share
        print(f"{sp:>13d}x {share:>11.0%} {saving:>23.1%}")
    print()
    print("  ==> The compute claim is explicitly scoped to one-at-a-time folding, which is how")
    print("      every notebook in this project ran ESMFold. If a reviewer batches, savings shrink")
    print("      roughly in proportion to the fold share. State the scope; do not claim it holds")
    print("      under arbitrary batching.")

    hr()
    print("BOTTOM LINE FOR THE PAPER")
    hr()
    if not sig:
        print("  §3k's savings numbers are methodologically sound on the axis that mattered most")
        print("  (fold cost does not depend on collapse status), and should be cited as-is, with")
        print("  the one-at-a-time-folding scope stated explicitly.")
    else:
        adj = "downward" if mean(ct) < mean(ht) else "upward"
        print(f"  §3k's savings need a {adj} adjustment of roughly {abs(1-ratio):.0%} on the folding")
        print(f"  component. Recompute before citing a headline percentage.")


if __name__ == "__main__":
    main()
