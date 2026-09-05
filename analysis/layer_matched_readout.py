"""
Consolidated read-out of the five matched-alpha_rel layer notebooks (49-53).

Three questions the raw summary CSVs do not answer on their own:

  1. Which models were actually RESPONSIVE? A layer comparison on an inert model compares two
     non-effects, and its null says nothing about layer depth.
  2. For the one model that moved a long way (Mistral-Prot), WHICH DIRECTION did it move, and is
     the pLDDT gain an artifact of degenerate repetition plus increased length?
  3. Does NB49 reproduce §1g's ProtGPT2 layer-12 result at the identical absolute push?

Zero-GPU, stdlib only (no pandas in the local env).
Run:  PYTHONIOENCODING=utf-8 python analysis/layer_matched_readout.py
"""

import collections
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = ["protgpt2", "zymctrl", "rita", "mistralprot", "progen2"]
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def load(model, kind):
    path = os.path.join(HERE, f"{model}_layer_matched_{kind}.csv")
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")


def hr(c="=", n=98):
    print(c * n)


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = mean(rx), mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


# ------------------------------------------------------------------ 1. cross-model overview
hr()
print("1. WHAT EACH MODEL DID  (matched alpha_rel; CONTROL vs its most extreme dose)")
hr()
print(f"{'model':<14}{'layers':<9}{'doses (a_rel)':<22}{'ctrl coll':>10}{'ctrl pLDDT':>11}"
      f"{'ctrl ent':>9}{'verdict':>14}")
print("-" * 98)

overview = {}
for m in MODELS:
    summ = load(m, "summary")
    cal = load(m, "calibration")[0]
    tests = load(m, "layer_tests")
    ctrl = next(r for r in summ if r["condition"] == "CONTROL")
    alphas = sorted({float(t["alpha_rel"]) for t in tests})
    any_layer_sig = any(t["significant"] == "True" for t in tests)
    overview[m] = dict(summ=summ, cal=cal, tests=tests, ctrl=ctrl,
                       any_layer_sig=any_layer_sig)
    print(f"{m:<14}{cal['layer_early']+'/'+cal['layer_late']:<9}"
          f"{', '.join(f'{a:.2f}' for a in alphas):<22}"
          f"{float(ctrl['collapse_rate'])*100:9.1f}%{float(ctrl['mean_plddt']):11.2f}"
          f"{float(ctrl['mean_entropy']):9.3f}"
          f"{('LAYER DIFF' if any_layer_sig else 'no layer diff'):>14}")

# ------------------------------------------------------------------ 2. alpha_rel machinery check
print()
hr()
print("2. IS THE alpha_rel MACHINERY SOUND? (this run's ||h|| vs NB34's audit)")
hr()
NB34 = {"protgpt2": (2921.53, 3751.99), "zymctrl": (51.14, 82.66),
        "rita": (23.69, 43.03), "mistralprot": (56439.23, 252032.61),
        "progen2": (None, None)}
print(f"{'model':<14}{'||h|| early':>14}{'NB34':>12}{'diff':>8}"
      f"{'||h|| late':>14}{'NB34':>12}{'diff':>8}{'old asym':>10}")
print("-" * 98)
for m in MODELS:
    cal = overview[m]["cal"]
    he, hl = float(cal["h_early"]), float(cal["h_late"])
    oe, ol = float(cal["old_alpha_rel_early"]), float(cal["old_alpha_rel_late"])
    r34e, r34l = NB34[m]
    de = f"{(he/r34e-1)*100:+.1f}%" if r34e else "  new"
    dl = f"{(hl/r34l-1)*100:+.1f}%" if r34l else "  new"
    print(f"{m:<14}{he:14.2f}{(r34e if r34e else 0):12.2f}{de:>8}"
          f"{hl:14.2f}{(r34l if r34l else 0):12.2f}{dl:>8}{oe/ol:9.2f}x")
print()
print("The 'old asym' column is what the ORIGINAL absolute-norm runs actually did: the early layer")
print("received that many times more relative push than the late layer. Every value > 1 means the")
print("early layer was over-pushed, in exactly the direction of the retracted early>late trend.")

# ------------------------------------------------------------------ 3. Mistral direction of effect
print()
hr()
print("3. MISTRAL-PROT: the effect is LARGE and points the WRONG WAY")
hr()
rows = load("mistralprot", "sequences")
by = collections.defaultdict(list)
for r in rows:
    by[r["condition"]].append(r)

print(f"{'condition':<12}{'collapse':>9}{'pLDDT':>8}{'entropy':>9}{'len':>7}"
      f"{'top-AA':>8}{'top-AA %':>10}{'K frac':>8}")
print("-" * 98)
order = ["CONTROL", "L2_D1x", "L2_D2x", "L2_D4x", "L7_D1x", "L7_D2x", "L7_D4x"]
mist_rows = []
for cond in order:
    g = by[cond]
    seqs = [r["sequence"] for r in g]
    lens = [len(s) for s in seqs]
    plddt = mean([float(r["plddt"]) for r in g])
    ent = mean([float(r["entropy"]) for r in g])
    coll = mean([int(r["collapse"]) for r in g])
    joined = "".join(seqs)
    cnt = collections.Counter(joined)
    top_aa, top_n = cnt.most_common(1)[0]
    kfrac = cnt["K"] / len(joined)
    mist_rows.append((cond, coll, plddt, ent, mean(lens), top_aa, top_n / len(joined), kfrac))
    print(f"{cond:<12}{coll*100:8.1f}%{plddt:8.2f}{ent:9.3f}{mean(lens):7.1f}"
          f"{top_aa:>8}{top_n/len(joined)*100:9.1f}%{kfrac*100:7.1f}%")

print()
print("Steering made the output MORE repetitive (entropy 3.50 -> 0.38) and roughly DOUBLED its")
print("length, and ESMFold rewarded it: collapse 94% -> 0%, pLDDT 42 -> 73. The L2_D4x output is")
print("~93% lysine. A poly-K tract is a textbook intrinsically-disordered, non-folding sequence;")
print("ESMFold confidently calls it well-structured.")

# length is a real confound here, so quantify it rather than hand-wave
allr = [r for r in rows]
sp_len = spearman([len(r["sequence"]) for r in allr], [float(r["plddt"]) for r in allr])
sp_ent = spearman([float(r["entropy"]) for r in allr], [float(r["plddt"]) for r in allr])
print()
print(f"Across all {len(allr)} Mistral sequences:  Spearman(length, pLDDT) = {sp_len:+.3f}")
print(f"                                          Spearman(entropy, pLDDT) = {sp_ent:+.3f}")
print("Both length and repetition track pLDDT here, and steering moved both at once -- so this")
print("run cannot separate 'poly-K folds confidently' from 'longer sequences score higher'.")
print("State it as a joint confound, not as a clean single-cause result.")

# ------------------------------------------------------------------ 4. ProtGPT2 replication
print()
hr()
print("4. DOES NB49 REPRODUCE §1g's ProtGPT2 RESULT AT THE IDENTICAL PUSH?")
hr()
g1g = {"CONTROL": (0.620, 57.84, 2.691), "L12_1x": (0.600, 58.14, 2.425),
       "L12_2x": (0.900, 42.66, 2.693), "L30_1x": (0.500, 59.80, 2.986),
       "L30_2x": (0.700, 56.30, 2.756)}
now = {r["condition"]: (float(r["collapse_rate"]), float(r["mean_plddt"]),
                        float(r["mean_entropy"])) for r in overview["protgpt2"]["summ"]}
pairs = [("CONTROL", "CONTROL"), ("L12_1x", "L12_D1x"), ("L12_2x", "L12_D2x"),
         ("L30_1x", "L30_D1x"), ("L30_2x", "L30_D2x")]
print(f"{'condition':<12}{'§1g coll':>10}{'NB49 coll':>11}{'§1g pLDDT':>11}{'NB49 pLDDT':>12}"
      f"{'d pLDDT':>10}{'verdict':>14}")
print("-" * 98)
for old_k, new_k in pairs:
    oc, op, oe = g1g[old_k]
    nc, np_, ne = now[new_k]
    d = np_ - op
    ok = abs(d) < 6.0 and abs(nc - oc) < 0.20
    print(f"{old_k:<12}{oc*100:9.1f}%{nc*100:10.1f}%{op:11.2f}{np_:12.2f}{d:+10.2f}"
          f"{('reproduces' if ok else 'DIVERGES'):>14}")

cal = overview["protgpt2"]["cal"]
print()
print(f"Push magnitude was identical by construction: matched_norm at L12 = "
      f"{float(cal['matched_norm_early']):.3f} vs §1g's 583.998.")
print(f"But the vector DIRECTION differs -- this run's raw ||v_L|| at L12 = "
      f"{float(cal['raw_v_norm_early']):.2f}, §1g's was 439.62.")
print()
print("Only L12 at 2x diverges; both L30 conditions and both 1x conditions reproduce. So this is")
print("NOT a broken pipeline -- it is one condition failing to replicate. Two candidate causes,")
print("which this run CANNOT separate:")
print("  (a) NB49 extracts the vector with a forward hook on block L, whereas 24/26/27/41/47 used")
print("      out.hidden_states[L]. In HF, hidden_states[L] is the output of block L-1, so every")
print("      earlier notebook built its vector ONE BLOCK EARLIER than it injected it. NB49 is")
print("      self-consistent; it is also therefore not a like-for-like repeat.")
print("  (b) §1g's 90% at L12_2x is a single N=50 point sitting on a steep dose-response cliff")
print("      (§1n: slope k=9.48), where a small change in vector direction can swing the outcome.")
print()
print("Decisive follow-up, ~40 min on a T4: build BOTH vectors from ONE pool, report their cosine")
print("similarity, and run 2x steering with each at N=50. That isolates (a) from (b).")

# ------------------------------------------------------------------ 5. bottom line
print()
hr()
print("5. WHAT THE LAYER QUESTION ACTUALLY GOT ANSWERED")
hr()
verdicts = {
    "protgpt2": ("INERT", "no condition moved; says nothing about layers"),
    "zymctrl": ("INERT", "3x did not reach significance here (§1h-THRESHOLD's 86% -> 76%)"),
    "rita": ("INERT", "consistent with §1k-REPAIRED; still nothing at alpha_rel 3.2"),
    "progen2": ("RESPONDED", "moved at 8x on the early layer, NO layer difference -> real null"),
    "mistralprot": ("RESPONDED", "moved hard at both layers, early >> late -> LAYER DIFFERENCE"),
}
for m in MODELS:
    v, note = verdicts[m]
    print(f"  {m:<14}{v:<12}{note}")
print()
print("Informative comparisons: 2 of 5 (ProGen2, Mistral-Prot). The other three were inert, which")
print("is a statement about dose, not about layer depth.")
print()
print("So the early>late claim does NOT come back in general -- but it is no longer unanimously")
print("dead either. Mistral-Prot shows a large, Holm-corrected early>late difference at genuinely")
print("matched relative push. Two caveats that must travel with it: the 'damage' is a pLDDT")
print("INCREASE via repetition collapse, and Mistral is an 8-expert MoE where a large constant")
print("injected at block 2 plausibly disrupts top-1 expert routing for every downstream block.")


if __name__ == "__main__":
    pass
