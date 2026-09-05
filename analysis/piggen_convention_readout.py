"""
NB58 read-out: does the extraction convention explain 1r-D, and does 1r-C survive?

Short answer: the convention does NOT explain it, 1r-C does NOT survive, and the run finds two
things nobody was looking for -- the steering vector's direction is only moderately reproducible,
and two NEARLY ORTHOGONAL vectors produce identical steering outcomes.

NOTE: scipy.stats is blocked on the local machine by an Application Control policy, so Fisher's
exact test is implemented here directly. Verified against the p-values NB58 itself produced with
real scipy (see the assertion at the bottom).

Run:  PYTHONIOENCODING=utf-8 python analysis/piggen_convention_readout.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _logcomb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher(a, b, c, d):
    # two-sided Fisher exact on [[a, b], [c, d]]
    n, r1, c1 = a + b + c + d, a + b, a + c
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)

    def pr(x):
        return math.exp(_logcomb(r1, x) + _logcomb(n - r1, c1 - x) - _logcomb(n, c1))

    p_obs = pr(a)
    return min(1.0, sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p_obs * (1 + 1e-9)))


def hr(c="=", n=94):
    print(c * n)


cos = load("piggen_vector_cosines.csv")
summ = {r["condition"]: r for r in load("piggen_convention_summary.csv")}
conv = load("piggen_convention_contrasts.csv")
lay = load("piggen_convention_layer_tests.csv")

hr()
print("1. THE CONVENTION GIVES A NEARLY ORTHOGONAL VECTOR -- and it changes nothing")
hr()
cv = {int(r["layer"]): float(r["value"]) for r in cos if r["metric"] == "convention"}
print(f"  cosine(v_hook, v_legacy)   layer 1: {cv[1]:+.4f}     layer 3: {cv[3]:+.4f}")
print(f"  ProtGPT2 (36 layers, 1q-C): +0.9985")
print()
print("The prediction was right about the geometry: on a 4-layer model, extracting one block")
print("earlier gives an almost unrelated direction (cosine 0.16-0.22, i.e. ~77-80 degrees apart),")
print("where on 36-layer ProtGPT2 the two were essentially identical.")
print()
print("But the steering outcome is the same either way:")
print(f"  {'layer':>6}{'hook':>18}{'legacy':>18}{'p_collapse':>12}{'p_pLDDT':>10}")
for r in conv:
    print(f"  {r['layer']:>6}{float(r['hook_rate'])*100:>10.0f}% /{float(r['hook_plddt']):>6.2f}"
          f"{float(r['legacy_rate'])*100:>10.0f}% /{float(r['legacy_plddt']):>6.2f}"
          f"{float(r['p_collapse']):>12.3g}{float(r['p_plddt']):>10.3g}")
print()
print("  ==> 1r-D is NOT explained by the extraction convention.")
print()
print("This is also, by accident, the cleanest test of 1l's 'magnitude dominates direction' the")
print("project has: two SEPARATELY CONSTRUCTED, nearly orthogonal repetition vectors, at matched")
print("alpha_rel, do the same thing. 1l needed artificial random/shuffled vectors to make that")
print("point; this is two plausible constructions of the same intended direction.")

hr()
print("2. THE STEERING VECTOR IS ONLY MODERATELY REPRODUCIBLE -- new, and it matters")
hr()
print(f"  {'convention':>12}{'layer':>7}{'split-half cosine':>20}{'sd':>8}{'worst split':>13}")
for r in cos:
    if r["metric"] == "convention":
        continue
    c = r["metric"].replace("split_half_", "")
    print(f"  {c:>12}{r['layer']:>7}{float(r['value']):>20.4f}"
          f"{float(r['sd']):>8.4f}{float(r['min']):>13.4f}")
sh = [float(r["value"]) for r in cos if r["metric"] != "convention"]
mins = [float(r["min"]) for r in cos if r["metric"] != "convention"]
print()
print(f"Two disjoint halves of the SAME pool give directions {min(sh):.2f}-{max(sh):.2f} cosine")
print(f"apart on average, and as low as {min(mins):.2f} on individual splits (~66 degrees).")
print()
print("The project scales every vector to a fixed alpha_rel, which controls MAGNITUDE. Nothing")
print("has ever controlled DIRECTION, and this is the first measurement of how much it wanders.")
print("It is the real explanation for 1r-D: each run draws a somewhat different vector.")
print()
print("Note the ordering that matters: the convention cosine (0.16-0.22) is FAR BELOW the")
print("split-half cosine (0.74-0.84). So the two conventions are not merely noisy versions of one")
print("another -- they are genuinely different directions, and they still behave identically.")

hr()
print("3. p-IgGen's LAYER-1 COLLAPSE RATE ACROSS FOUR RUNS, all at alpha_rel 0.400")
hr()
runs = [("41  (legacy)", 0, 50), ("57  (hook)", 17, 100),
        ("58  (hook)", 9, 100), ("58  (legacy)", 9, 100)]
print(f"  {'run':>14}{'collapsed':>12}{'rate':>9}")
for name, k, n in runs:
    print(f"  {name:>14}{k:>7}/{n:<4}{k/n:>8.0%}")
print()
print("  pairwise Fisher exact:")
for i in range(len(runs)):
    for j in range(i + 1, len(runs)):
        (n1, k1, t1), (n2, k2, t2) = runs[i], runs[j]
        p = fisher(k1, t1 - k1, k2, t2 - k2)
        flag = "  <-- differ" if p < 0.05 else ""
        print(f"    {n1:>14} vs {n2:<14} p = {p:.4f}{flag}")
print()
print("Only 41 stands apart, and only against the 17% run. The 0%/9%/9%/17% spread is what a")
print("moderately unstable vector plus a small true effect looks like.")

hr()
print("4. DOES 1r-C SURVIVE? NO.")
hr()
print(f"  {'convention':>12}{'early':>16}{'late':>16}{'d pLDDT':>10}{'p_pLDDT':>11}{'verdict':>10}")
for r in lay:
    print(f"  {r['convention']:>12}{float(r['early_rate'])*100:>8.0f}% /{float(r['early_plddt']):>6.2f}"
          f"{float(r['late_rate'])*100:>8.0f}% /{float(r['late_plddt']):>6.2f}"
          f"{float(r['diff_plddt']):>10.2f}{float(r['p_plddt']):>11.3g}"
          f"{('SIG' if r['layer_difference'] == 'True' else 'ns'):>10}")
print()
print("  57 (hook, N=100):  +15 pts collapse, -4.96 pLDDT, p = 2.0e-04   SIG")
print("  58 (hook, N=100):   +5 pts collapse, -1.68 pLDDT, p = 7.3e-02   ns")
print("  58 (legacy,N=100):  +8 pts collapse, -3.64 pLDDT, p = 4.8e-04   SIG")
print()
print("Same model, same dose, same length, same N, and the layer difference is significant in two")
print("runs out of three -- with the two significant ones using DIFFERENT conventions and the two")
print("hook runs disagreeing with each other. That is not a robust finding.")
print()
print("  ==> 1r-C must be RETRACTED. p-IgGen does not show a reliable early>late layer difference.")
print("      The pre-registered reading for 'significant under only one convention' was 'a serious")
print("      caveat'; the actual pattern is worse than that -- it is not stable across replications")
print("      of the SAME convention.")

hr()
print("5. WHAT SURVIVES")
hr()
print("  * Mistral-Prot's layer difference stands (p_pLDDT 2.1e-16, effect +24.6 pLDDT). It is")
print("    two orders of magnitude larger than anything here and is not vulnerable to this.")
print("  * ProtGPT2 (1r-A) and ZymCTRL (1r-B): no layer difference. Unchanged.")
print("  * p-IgGen: RETRACTED -- no reliable layer difference.")
print("  * So the honest cross-model claim becomes: injection depth matters on ONE model of the")
print("    five testable (Mistral-Prot), and that model is an MoE where a large constant at block 2")
print("    plausibly disrupts expert routing. That is a much weaker claim than 1r suggested, and")
print("    it is the one the data supports.")
print()
print("  NEW, and worth more than the layer question:")
print("  * Steering-vector direction is only moderately reproducible (split-half cosine 0.74-0.84).")
print("  * Two nearly orthogonal vectors (cosine 0.16) steer identically -- independent support for")
print("    1l, on a second model, without artificial random directions.")

# Guard: the Fisher implementation here must agree with the scipy values NB58 printed.
_p = fisher(9, 91, 4, 96)
assert abs(_p - float(next(r["p_collapse"] for r in lay if r["convention"] == "hook"))) < 1e-6, \
    "local Fisher disagrees with the notebook's scipy value -- do not trust these p-values"
print()
print("(local Fisher implementation verified against NB58's own scipy output)")
