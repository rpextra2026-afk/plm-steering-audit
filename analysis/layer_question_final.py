"""
The early-vs-late layer question, FINAL: every model tested at a dose where it demonstrably
responds, at fixed generation length, with adequate N.

Supersedes layer_question_summary.py, which covered only the underpowered / free-length runs.

Sources (all in analysis/ except NB54, whose CSVs were never downloaded -- its numbers are quoted
from the notebook's saved output cells, see locked-results.md 1q):
  NB55  Mistral-Prot  fixed length, alpha_rel 0.80   mistralprot_fixedlen_*.csv
  NB56  ProtGPT2      fixed length, N=150/arm        protgpt2_layer_n150_*.csv
  NB57  ZymCTRL       screen -> confirm, N=100/arm   zymctrl_*.csv
  NB57  p-IgGen       screen -> confirm, N=100/arm   piggen_*.csv

Zero-GPU, stdlib only.
Run:  PYTHONIOENCODING=utf-8 python analysis/layer_question_final.py
"""

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        return None
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def hr(c="=", n=100):
    print(c * n)


# model, layers, dose alpha_rel, N/arm, early-vs-control pLDDT + p, late-vs-control pLDDT + p,
# layer diff pLDDT, p, verdict
FINAL = [
    ("ProtGPT2", "12/30", 36, 0.200, 150, -4.42, 0.000264, -1.73, 0.2327, -2.69, 0.005238,
     "no (below floor)"),
    ("ZymCTRL", "12/30", 36, 0.400, 100, -3.64, 0.01222, -2.55, 0.04739, -1.09, 0.7058, "no"),
    ("p-IgGen", "1/3", 4, 0.400, 100, -6.00, 4.728e-07, -1.03, 0.09515, -4.96, 0.0002031,
     "YES"),
    ("Mistral-Prot", "2/7", 8, 0.800, 50, +26.05, 2.07e-17, +1.50, 0.992, +24.55, 2.14e-16,
     "YES"),
]

hr()
print("THE LAYER QUESTION, ANSWERED PROPERLY")
print("every model at a dose where it responds, fixed length, adequate N")
hr()
print("{:14s}{:>7}{:>7}{:>7}{:>6}{:>20}{:>20}{:>20}".format(
    "model", "layers", "depth", "a_rel", "N", "early vs ctrl", "late vs ctrl", "early vs late"))
print("-" * 100)
for (m, lay, dep, a, n, de, pe, dl, pl, dd, pd_, verdict) in FINAL:
    print("{:14s}{:>7}{:>7}{:7.2f}{:6d}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}".format(
        m, lay, dep, a, n,
        f"{de:+.2f}", f"p={pe:.2g}", f"{dl:+.2f}", f"p={pl:.2g}",
        f"{dd:+.2f}", f"p={pd_:.2g}"))
print()
print("{:14s}{:>20}".format("model", "LAYER DIFFERENCE?"))
print("-" * 40)
for (m, lay, dep, a, n, de, pe, dl, pl, dd, pd_, verdict) in FINAL:
    print("{:14s}{:>20}".format(m, verdict))
print("  RITA          untestable -- flat to alpha_rel 3.2 (16x ProtGPT2's threshold)")
print("  ProGen2       untestable -- inert, and 22-residue output where pLDDT means least")

hr()
print("1. ProtGPT2 -- the N=50 gap did NOT replicate")
hr()
print("           NB54 (N=50)   NB56 (N=150)")
print("  collapse gap   +18.0 pts      +3.3 pts   CI [-6.0, +12.7]  p=0.59")
print("  pLDDT gap       -4.80         -2.69      CI [-5.09, -0.21] p=0.0052")
print()
print("Textbook winner's curse: the N=50 estimate was inflated. The collapse gap is gone.")
print("The pLDDT gap survives Holm (p=0.005) but is -2.69, BELOW the pre-registered 3-point")
print("materiality floor -- so the pre-registered verdict is NOT SUPPORTED, and that is what")
print("gets reported. This is exactly what pre-registration is for: without the floor written")
print("down in advance, p=0.005 would have been claimed as a win.")
print()
print("A tempting misreading to avoid: layer 12 moved vs control (pLDDT -4.42, p=0.00026) while")
print("layer 30 did not (-1.73, p=0.23). That is NOT evidence of a layer difference -- the")
print("difference between 'significant' and 'not significant' is not itself significant. The")
print("direct contrast is the pre-registered test, and it gives -2.69 [-5.09, -0.21].")

hr()
print("2. ZymCTRL -- responsive at last, and no layer difference")
hr()
print("The screen found alpha_rel 0.40 responsive; at N=100 the early layer moved (collapse")
print("+15 pts p=0.026, pLDDT -3.64 p=0.012). So the comparison was LIVE, and the layer null")
print("is interpretable: +7.0 pts collapse (p=0.32), -1.09 pLDDT (p=0.71).")
print()
print("Caveat on the screen: at N=25/rung it was non-monotonic -- 2x flagged responsive but")
print("3x, 4x and 6x did not. That is screen noise. The pick was nevertheless validated by the")
print("confirm run's own control comparison, which is what the two-phase design is for.")

hr()
print("3. p-IgGen -- a clean, properly-controlled LAYER DIFFERENCE")
hr()
print("  CONTROL   1% collapse / 74.45 pLDDT")
print("  layer 1  17% collapse / 68.46 pLDDT   vs control: p=7.5e-05 / 4.7e-07")
print("  layer 3   2% collapse / 73.42 pLDDT   vs control: p=1.0    / 0.095")
print("  early vs late: +15 pts collapse p=0.00041, -4.96 pLDDT p=0.00020, BOTH Holm-significant")
print()
print("The early layer damages; the late layer does essentially nothing. Dose chosen on")
print("independent data (screen, seeds 7000+), tested on fresh samples (seeds 9000+).")
print()
print("Note the screen's value here: 4x and 8x drove collapse to 88% and 100% -- saturated, where")
print("a layer comparison would be uninformative. Taking the LOWEST responsive dose avoided that.")

hr()
print("4. THE UNRESOLVED PROBLEM -- p-IgGen contradicts NB41")
hr()
print("NB41 tested p-IgGen at the SAME alpha_rel (0.40), at effectively the same fixed length")
print("(every sequence exactly 49.0 residues), and found 0/50 collapse at layer 1.")
print("NB57 finds 17/100. If the true rate were 0.17, P(0 of 50) is about 1e-4 -- so this is")
print("NOT sampling noise. Something about the intervention differs.")
print()
print("Two candidates, and this data cannot separate them:")
print("  (a) THE EXTRACTION CONVENTION. 1q-C established the hidden_states[L] vs hook-on-block-L")
print("      off-by-one is harmless -- but that was measured on ProtGPT2, which has 36 layers,")
print("      where one block is 2.8% of the network. p-IgGen has FOUR. One block is 25%.")
print("      The cosine similarity between the two vectors has never been measured on p-IgGen.")
print("      >>> 1q-C's 'harmless' conclusion is a ProtGPT2 result and must NOT be generalised")
print("      >>> to shallow models until this is checked.")
print("  (b) A different candidate pool (200 in NB41 vs 150 in NB57) giving a different vector.")
print()
print("Cheap decisive test, minutes and no folding: build both vectors on p-IgGen from one pool")
print("and print cosine(v_hook, v_legacy) at layers 1 and 3.")

hr()
print("5. THE PATTERN WORTH NOTICING (hypothesis, not a claim)")
hr()
print("  layer difference FOUND:     p-IgGen (4 layers), Mistral-Prot (8 layers)")
print("  layer difference NOT found: ProtGPT2 (36 layers), ZymCTRL (36 layers)")
print()
print("The two models where injection depth matters are the two SHALLOW ones. Relative depth is")
print("similar across all four (~25-33%), so this is about total depth, not position: in a")
print("36-layer model, layers 12 and 30 are both 'middle' and functionally similar, and 24 blocks")
print("remain to absorb the perturbation. In a 4-layer model, block 1 and block 3 do genuinely")
print("different jobs and only 3 blocks follow.")
print()
print("Speculative. Two models per side is not a result. But it is a coherent reading of why the")
print("original early>late claim looked real on p-IgGen and nowhere else -- and it predicts the")
print("effect should track total depth, which is testable.")

hr()
print("6. WHAT TO WRITE")
hr()
print("The retracted early>late claim does NOT come back as a general finding. It comes back")
print("as a MODEL-SPECIFIC one, on exactly the model that originally carried it (p-IgGen), plus")
print("Mistral-Prot -- and it stays dead on the two deep models.")
print()
print("Do not write 'un-retracted' until the NB41 contradiction in section 4 is resolved. Until")
print("then the honest sentence is:")
print()
print("  'At matched relative push and fixed generation length, injection depth affects")
print("   structural damage in the two shallow models tested (p-IgGen, Mistral-Prot) and not in")
print("   the two deep ones (ProtGPT2, ZymCTRL). The p-IgGen result conflicts with an earlier")
print("   run of ours at the same dose, which we have not resolved.'")
