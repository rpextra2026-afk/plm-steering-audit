"""
Pooled statistical re-analysis of already-locked collapse counts.

No new experiments and no GPU: every count below is transcribed from a saved notebook output
cell via notes/locked-results.md. The point is that several of this project's findings are
"consistent in direction but individually non-significant at N=50", which is a weak way to
report a result that standard meta-analytic pooling can state properly.

Two analyses:

  A. Early-vs-late layer sensitivity across the four informative models (locked-results SS1k).
     Reported there as "4-for-4 in direction, only p-IgGen individually significant".
     Pooled here with Mantel-Haenszel, plus an exact sign test.

  B. Decoding-only interventions vs CONTROL (locked-results SS1c). Reported there as
     "a consistent but not individually significant trend" across three conditions, each of
     whose CI overlaps CONTROL's. Pooled here as a single "any decoding intervention" arm.

IMPORTANT CAVEAT on analysis A, which must travel with any number it produces: pooling assumes
the strata are measuring the same underlying effect. Per notes/locked-results.md SS1k-FLAG, the
early and late layers were NOT pushed at matched relative strength -- the early layer received a
1.44x-3.23x larger push relative to its own residual-stream scale in every model. So a
significant pooled result here does NOT establish that early layers are intrinsically more
fragile; it establishes that the early-layer CONDITION produced more collapse, which is at
present partly a statement about push size. Notebook 34 measures the relevant norms and settles
which reading is correct. Until then, treat A as "the pooled version of a confounded
comparison", not as confirmation.

Run:  python analysis/meta_analysis.py
"""

import math
from itertools import combinations

# --------------------------------------------------------------------------------------
# Data. (collapsed, total) per condition, transcribed from notes/locked-results.md.
# --------------------------------------------------------------------------------------

# A. Early vs late layer, both at nominal 1x reference strength (SS1g/SS1h/SS1i/SS1k).
LAYER_STRATA = [
    # model,               early (collapsed, n),  late (collapsed, n),  push asymmetry (SS1k-FLAG)
    ("ProtGPT2  L12/L30",  (30, 50),              (25, 50),             1.78),
    ("ZymCTRL   L12/L30",  (43, 50),              (36, 50),             1.77),
    ("p-IgGen   L1/L3",    (50, 50),              (16, 50),             3.23),
    ("RITA      L3/L11",   (48, 50),              (47, 50),             1.80),
]

# B. Decoding-only interventions vs a shared CONTROL (SS1c).
DECODING_CONTROL = (22, 50)
DECODING_ARMS = [
    ("REP_PENALTY_1.2",   (29, 50)),
    ("REP_PENALTY_1.5",   (35, 50)),
    ("NO_REPEAT_NGRAM_3", (30, 50)),
]


# --------------------------------------------------------------------------------------
# Minimal stats, implemented directly so this script has no dependencies beyond stdlib.
# --------------------------------------------------------------------------------------

def norm_sf(z):
    """Upper-tail probability of the standard normal."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def two_sided_z_p(z):
    return 2.0 * norm_sf(abs(z))


def log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_two_sided(a, b, c, d):
    """Two-sided Fisher exact test on [[a,b],[c,d]], by summing all tables at most as probable."""
    row1, row2 = a + b, c + d
    col1 = a + c
    total = row1 + row2
    log_denom = log_comb(total, col1)

    def logp(x):
        return log_comb(row1, x) + log_comb(row2, col1 - x) - log_denom

    lo = max(0, col1 - row2)
    hi = min(col1, row1)
    p_obs = logp(a)
    tol = 1e-9
    total_p = 0.0
    for x in range(lo, hi + 1):
        lp = logp(x)
        if lp <= p_obs + tol:
            total_p += math.exp(lp)
    return min(1.0, total_p)


def binom_test_two_sided(k, n, p=0.5):
    """Exact binomial test by summing outcomes no more probable than the observed one."""
    p_obs = math.exp(log_comb(n, k) + k * math.log(p) + (n - k) * math.log(1 - p))
    tol = 1e-9
    total = 0.0
    for i in range(n + 1):
        pi = math.exp(log_comb(n, i) + i * math.log(p) + (n - i) * math.log(1 - p))
        if pi <= p_obs + tol:
            total += pi
    return min(1.0, total)


def wilson_ci(k, n, z=1.959963985):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def mantel_haenszel(strata, correction=0.0):
    """
    Mantel-Haenszel pooled odds ratio + the MH chi-square test.
    Each stratum: ((a, n1), (c, n2)) where a = events in arm 1, c = events in arm 2.
    `correction` adds a constant to every cell (Haldane-Anscombe) to tolerate zero cells.
    """
    num = den = 0.0
    s_num = s_den = 0.0
    sum_a = sum_ea = sum_va = 0.0

    for (a, n1), (c, n2) in strata:
        a = a + correction
        b = (n1 - a) + 2 * correction if correction else n1 - a
        c = c + correction
        d = (n2 - c) + 2 * correction if correction else n2 - c
        # keep row totals consistent after correction
        b = max(b, 1e-12)
        d = max(d, 1e-12)
        t = a + b + c + d

        num += a * d / t
        den += b * c / t

        # MH chi-square components, computed on the corrected table
        r1, r2 = a + b, c + d
        c1 = a + c
        sum_a += a
        sum_ea += r1 * c1 / t
        if t > 1:
            sum_va += (r1 * r2 * c1 * (t - c1)) / (t * t * (t - 1))

    if den == 0 or num == 0:
        return None

    or_mh = num / den

    # Robins-Breslow-Greenland variance of log(OR_MH)
    p_sum = q_sum = pr_sum = ps_qr_sum = qs_sum = 0.0
    for (a, n1), (c, n2) in strata:
        a = a + correction
        b = (n1 - a) + 2 * correction if correction else n1 - a
        c = c + correction
        d = (n2 - c) + 2 * correction if correction else n2 - c
        b = max(b, 1e-12)
        d = max(d, 1e-12)
        t = a + b + c + d
        P = (a + d) / t
        Q = (b + c) / t
        R = a * d / t
        S = b * c / t
        p_sum += P * R
        q_sum += Q * S
        pr_sum += P * S + Q * R
        ps_qr_sum += R
        qs_sum += S

    var_log_or = (p_sum / (2 * ps_qr_sum ** 2)
                  + pr_sum / (2 * ps_qr_sum * qs_sum)
                  + q_sum / (2 * qs_sum ** 2))
    se = math.sqrt(var_log_or)
    lo = math.exp(math.log(or_mh) - 1.959963985 * se)
    hi = math.exp(math.log(or_mh) + 1.959963985 * se)

    chi2 = ((abs(sum_a - sum_ea) - 0.5) ** 2) / sum_va if sum_va > 0 else 0.0
    p = two_sided_z_p(math.sqrt(chi2))
    return {"or": or_mh, "ci": (lo, hi), "chi2": chi2, "p": p, "se_log_or": se}


def hr(char="=", n=88):
    print(char * n)


# --------------------------------------------------------------------------------------
# Analysis A -- early vs late layer sensitivity
# --------------------------------------------------------------------------------------

def analysis_a():
    hr()
    print("A. EARLY vs LATE LAYER SENSITIVITY AT NOMINAL 1x  (locked-results SS1g-SS1k)")
    hr()
    print("Currently reported as: '4-for-4 in direction, only individually significant in p-IgGen'.")
    print()
    print(f"{'Model':22s} {'early':>12s} {'late':>12s} {'diff':>8s} {'OR':>8s} "
          f"{'Fisher p':>10s} {'push asym':>10s}")
    print("-" * 88)

    strata = []
    n_correct_direction = 0
    for name, (ce, ne), (cl, nl) in [(s[0], s[1], s[2]) for s in LAYER_STRATA]:
        pass  # placeholder, real loop below

    for name, early, late, asym in LAYER_STRATA:
        ce, ne = early
        cl, nl = late
        pe, pl = ce / ne, cl / nl
        diff = pe - pl
        if diff > 0:
            n_correct_direction += 1
        a, b = ce, ne - ce
        c, d = cl, nl - cl
        p_fisher = fisher_exact_two_sided(a, b, c, d)
        if b == 0 or c == 0:
            or_str = "inf"
        else:
            or_str = f"{(a * d) / (b * c):8.2f}"
        print(f"{name:22s} {ce:4d}/{ne:<3d}{pe:6.0%} {cl:4d}/{nl:<3d}{pl:6.0%} "
              f"{diff:+7.1%} {or_str:>8s} {p_fisher:10.4f} {asym:9.2f}x")
        strata.append(((ce, ne), (cl, nl)))

    print()
    print(f"Direction consistent in {n_correct_direction}/{len(LAYER_STRATA)} models.")

    sign_p = binom_test_two_sided(n_correct_direction, len(LAYER_STRATA), 0.5)
    print(f"Exact sign test (is 'always the same direction' more than chance?): p = {sign_p:.4f}")
    print("  -> with only 4 strata the sign test alone can never reach p<0.05 even at 4/4")
    print("     (the floor is p=0.125 two-sided), which is exactly why pooling the counts")
    print("     rather than the directions is the right move here.")
    print()

    mh = mantel_haenszel(strata, correction=0.0)
    mh_corr = mantel_haenszel(strata, correction=0.5)

    print("Mantel-Haenszel pooled odds ratio (early vs late, OR>1 = early collapses more):")
    if mh:
        print(f"  all 4 models, uncorrected : OR = {mh['or']:.2f}  "
              f"95% CI [{mh['ci'][0]:.2f}, {mh['ci'][1]:.2f}]   p = {mh['p']:.3g}")
    if mh_corr:
        print(f"  all 4 models, +0.5 corrected: OR = {mh_corr['or']:.2f}  "
              f"95% CI [{mh_corr['ci'][0]:.2f}, {mh_corr['ci'][1]:.2f}]   p = {mh_corr['p']:.3g}")
    print("  (p-IgGen contributes a zero cell -- 50/50 collapsed -- so the corrected row is the")
    print("   one to quote; the uncorrected row is shown only to confirm the correction is not")
    print("   doing the work.)")
    print()

    # Leave-one-out: does the whole result rest on p-IgGen?
    print("Leave-one-out sensitivity (does the pooled result depend on any single model?):")
    for i, (name, _, _, _) in enumerate(LAYER_STRATA):
        reduced = [s for j, s in enumerate(strata) if j != i]
        m = mantel_haenszel(reduced, correction=0.5)
        if m:
            print(f"  without {name:22s}: OR = {m['or']:5.2f}  "
                  f"95% CI [{m['ci'][0]:.2f}, {m['ci'][1]:.2f}]  p = {m['p']:.3g}")
    print()
    print("READ THIS WITH SS1k-FLAG. The 'push asym' column is the factor by which the early")
    print("layer was pushed harder relative to its own residual-stream scale. It is >1 in every")
    print("row, and largest exactly where the effect is largest (p-IgGen). A pooled OR computed")
    print("across strata that were not equally dosed measures dose as well as depth. Notebook 34")
    print("resolves which. Do NOT put this pooled number in the paper as evidence for intrinsic")
    print("early-layer fragility until it does.")
    print()


# --------------------------------------------------------------------------------------
# Analysis B -- decoding-only interventions
# --------------------------------------------------------------------------------------

def analysis_b():
    hr()
    print("B. DECODING-ONLY INTERVENTIONS vs CONTROL  (locked-results SS1c)")
    hr()
    print("Currently reported as: 'a consistent but not individually significant trend' --")
    print("all three CIs overlap CONTROL's at N=50.")
    print()

    cc, cn = DECODING_CONTROL
    clo, chi = wilson_ci(cc, cn)
    print(f"{'Condition':22s} {'collapsed':>12s} {'rate':>8s} {'95% Wilson CI':>20s} {'Fisher p vs CONTROL':>21s}")
    print("-" * 88)
    print(f"{'CONTROL':22s} {cc:6d}/{cn:<5d} {cc/cn:7.1%} {f'[{clo:.1%}, {chi:.1%}]':>20s} {'--':>21s}")

    tot_c = tot_n = 0
    for name, (k, n) in DECODING_ARMS:
        lo, hi = wilson_ci(k, n)
        p = fisher_exact_two_sided(k, n - k, cc, cn - cc)
        print(f"{name:22s} {k:6d}/{n:<5d} {k/n:7.1%} {f'[{lo:.1%}, {hi:.1%}]':>20s} {p:21.4f}")
        tot_c += k
        tot_n += n

    print("-" * 88)
    lo, hi = wilson_ci(tot_c, tot_n)
    p_pooled = fisher_exact_two_sided(tot_c, tot_n - tot_c, cc, cn - cc)
    print(f"{'ALL THREE POOLED':22s} {tot_c:6d}/{tot_n:<5d} {tot_c/tot_n:7.1%} "
          f"{f'[{lo:.1%}, {hi:.1%}]':>20s} {p_pooled:21.4f}")

    a, b, c, d = tot_c, tot_n - tot_c, cc, cn - cc
    or_pooled = (a * d) / (b * c) if b and c else float("inf")
    print()
    print(f"Pooled 'any decoding-level intervention' vs CONTROL: OR = {or_pooled:.2f}, "
          f"Fisher p = {p_pooled:.4f}")
    print()
    if p_pooled < 0.05:
        print("  ==> The pooled comparison DOES clear significance, unlike any individual arm.")
        print("      This upgrades SS1c from 'a directional trend we hedge' to a real, reportable")
        print("      finding: decoding-level repetition controls significantly WORSEN structural")
        print("      collapse relative to no intervention. That strengthens the control's actual")
        print("      job in the paper -- showing the steering effect is not just 'any change to")
        print("      the repetition profile' -- and it costs one sentence, no new GPU time.")
    else:
        print("  ==> Still not significant even pooled. Keep the existing hedge in SS1c; the")
        print("      honest statement is that the control cannot distinguish these from CONTROL")
        print("      at this N.")
    print()
    print("Caveat to state if this is used: the three arms share one CONTROL group, so they are")
    print("not independent strata and this is a pooled-arm comparison, not a meta-analysis.")
    print("It answers 'does any decoding-level intervention differ from control', which is the")
    print("question SS1c actually needs answered.")
    print()


if __name__ == "__main__":
    print()
    analysis_a()
    analysis_b()
    hr()
    print("All inputs transcribed from notes/locked-results.md. No new experiments were run.")
    hr()
