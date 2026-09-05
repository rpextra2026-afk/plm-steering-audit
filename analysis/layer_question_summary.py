"""
The early-vs-late layer question, every test the project has run at matched alpha_rel, in one
table -- plus a power analysis asking whether the nulls are "no effect" or "not enough N".

Sources:
  NB41  p-IgGen     free length   (piggen_layer_matched_summary.csv)
  NB49-53            free length   ({model}_layer_matched_layer_tests.csv)
  NB54  ProtGPT2    FIXED length  (quoted from the notebook's saved output cells; CSVs not
                                   downloaded from Kaggle -- see locked-results.md 1q)
  NB55  Mistral/ProGen2 FIXED len ({model}_fixedlen_layer_tests.csv)

Zero-GPU, stdlib + scipy.
Run:  PYTHONIOENCODING=utf-8 python analysis/layer_question_summary.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        return None
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def hr(c="=", n=104):
    print(c * n)


# ---------------------------------------------------------------- the consolidated table
# (model, run, alpha_rel, early%, late%, early pLDDT, late pLDDT, p_collapse, p_pLDDT,
#  responded?, layer difference?)
ROWS = []

# NB41 -- p-IgGen, free length. Both layers 0/50 at both doses; Fisher p = 1.0 both.
for a in (0.20, 0.40):
    ROWS.append(("p-IgGen", "41 free", a, 0.0, 0.0, None, None, 1.0, None, False, False))

# NB49-53, free length
for tag, model in [("protgpt2", "ProtGPT2"), ("zymctrl", "ZymCTRL"), ("rita", "RITA"),
                   ("mistralprot", "Mistral-Prot"), ("progen2", "ProGen2")]:
    rows = load(f"{tag}_layer_matched_layer_tests.csv")
    if not rows:
        continue
    # responsiveness verdicts as recorded in locked-results 1p
    responded = {"protgpt2": False, "zymctrl": False, "rita": False,
                 "mistralprot": True, "progen2": True}[tag]
    for r in rows:
        ROWS.append((model, "49-53 free", float(r["alpha_rel"]),
                     float(r["early_rate"]), float(r["late_rate"]),
                     float(r["early_plddt"]), float(r["late_plddt"]),
                     float(r["p_fisher_collapse"]), float(r["p_mwu_plddt"]),
                     responded, r["significant"] == "True"))

# NB54 -- ProtGPT2 at FIXED length, from the notebook's saved output
ROWS.append(("ProtGPT2", "54 FIXED", 0.200, 0.84, 0.66, 49.81, 54.61, 0.0634, 0.0281, True, False))
ROWS.append(("ProtGPT2", "54 FIXED", 0.400, 0.76, 0.70, 53.10, 53.86, 0.653, 0.749, True, False))

# NB55 -- Mistral / ProGen2 at FIXED length
for tag, model, responded in [("mistralprot", "Mistral-Prot", True), ("progen2", "ProGen2", False)]:
    rows = load(f"{tag}_fixedlen_layer_tests.csv")
    if not rows:
        continue
    for r in rows:
        ROWS.append((model, "55 FIXED", float(r["alpha_rel"]),
                     float(r["early_rate"]), float(r["late_rate"]),
                     float(r["early_plddt"]), float(r["late_plddt"]),
                     float(r["p_collapse"]), float(r["p_plddt"]),
                     responded, r["significant"] == "True"))

hr()
print("EVERY EARLY-vs-LATE TEST AT MATCHED alpha_rel")
hr()
print("{:14s}{:12s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>11s}{:>11s}{:>7s}{:>7s}".format(
    "model", "run", "a_rel", "early%", "late%", "eLDDT", "lLDDT", "p_coll", "p_pLDDT",
    "resp", "SIG"))
print("-" * 104)
last = None
for (m, run, a, er, lr, ep, lp, pc, pp, resp, sig) in ROWS:
    if last and last != (m, run):
        print()
    last = (m, run)
    epl = f"{ep:8.1f}" if ep is not None else "       -"
    lpl = f"{lp:8.1f}" if lp is not None else "       -"
    ppl = f"{pp:11.3g}" if pp is not None else "          -"
    print("{:14s}{:12s}{:8.2f}{:7.0f}%{:7.0f}%{:s}{:s}{:11.3g}{:s}{:>7s}{:>7s}".format(
        m, run, a, er * 100, lr * 100, epl, lpl, pc, ppl,
        "yes" if resp else "no", "YES" if sig else "-"))

print()
print("resp = did the model respond to steering AT ALL at that dose (vs its own control)?")
print("A layer comparison on a model that did not respond is comparing two non-effects: those")
print("rows are not evidence about layer depth, they are evidence about dose.")

# ---------------------------------------------------------------- power
hr()
print("WAS THE ProtGPT2 LAYER NULL 'NO EFFECT', OR 'NOT ENOUGH N'?")
hr()


def phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def phi_inv(p):
    # Acklam's rational approximation, plenty accurate here
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > ph:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def cohen_h(p1, p2):
    return abs(2 * math.asin(math.sqrt(p1)) - 2 * math.asin(math.sqrt(p2)))


def power(h, n, alpha):
    return phi(h * math.sqrt(n / 2.0) - phi_inv(1 - alpha / 2))


def n_for(h, alpha, target=0.80):
    return 2 * ((phi_inv(1 - alpha / 2) + phi_inv(target)) / h) ** 2


# The observed ProtGPT2 fixed-length layer gap at 1x
p_early, p_late, n_used = 0.84, 0.66, 50
h = cohen_h(p_early, p_late)
print(f"Observed at 1x, fixed length: L12 {p_early:.0%} collapse vs L30 {p_late:.0%}"
      f"  (gap {(p_early-p_late)*100:.0f} points)")
print(f"Cohen's h = {h:.4f}")
print()
print(f"{'alpha used':<28}{'power at n=50':>16}{'n/arm for 80% power':>24}")
print("-" * 70)
for label, a in [("0.05, uncorrected", 0.05), ("0.0250 (Holm, 2 tests)", 0.025),
                 ("0.0125 (Holm, 4 tests)", 0.0125)]:
    print(f"{label:<28}{power(h, n_used, a)*100:15.0f}%{n_for(h, a):24.0f}")

print()
print("So at N=50 with the correction actually applied, the chance of detecting a gap of exactly")
print("the size we observed was roughly ONE IN THREE. The null is not evidence of no effect --")
print("it is what an underpowered test looks like. Note the pLDDT readout agrees: MWU p = 0.0281")
print("at 1x, which is significant uncorrected and fails only under multiplicity correction.")

print()
print(f"{'design':<44}{'power (Holm, 2 tests)':>24}")
print("-" * 70)
for n in (50, 100, 150, 200):
    print(f"{'N=' + str(n) + ' per arm, single dose, 2 layers':<44}{power(h, n, 0.025)*100:23.0f}%")

hr()
print("WHERE THE OTHER MODELS ACTUALLY STAND")
hr()
print("The layer question can only be asked of a model that responds. Ranked by how strongly each")
print("responded at its best dose, at matched alpha_rel:")
print()
print("  Mistral-Prot  +26.1 pLDDT / -96 pts collapse   -> LAYER DIFFERENCE, p_pLDDT 2.1e-16")
print("  ProtGPT2       -7.1 pLDDT / +28 pts collapse   -> trend, p_pLDDT 0.028, fails correction")
print("  ProGen2        -3.4 pLDDT (fails correction)   -> cannot ask the question")
print("  ZymCTRL        -3.3 pLDDT (fails correction)   -> cannot ask the question")
print("  RITA           flat to alpha_rel 3.2           -> cannot ask the question")
print()
print("The ordering is the point: the size of the layer difference tracks the size of the overall")
print("response. That is exactly the pattern you would see if layer effects are real but need a")
print("large enough response to be detectable at N=50 -- and exactly what you would ALSO see if")
print("only Mistral has one. This data cannot separate those two readings.")
print()
print("Per-model reading of each null:")
print("  ZymCTRL  -- plausibly dose. It was trending (76% vs 62% control at alpha_rel 0.60,")
print("              p=0.19) and 1h-THRESHOLD once measured 86% there. Untested above 0.60.")
print("  RITA     -- probably genuine robustness. Flat to alpha_rel 3.2, which is 16x ProtGPT2's")
print("              own threshold. Not a dose problem.")
print("  ProGen2  -- unclear, and its natural output is 22 residues, where pLDDT means least.")
print("  p-IgGen  -- 0% collapse at both layers, both doses. Untested above alpha_rel 0.40.")
