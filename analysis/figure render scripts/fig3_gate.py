# fig3_gate.py -- Figure 3: "the gate is not the truth". THE WELD.
#
# This is the only figure that puts beats 3 and 7 in one frame. They are the
# same claim measured twice, at the two links of one chain:
#
#   (a)  the metric the published steering method OPTIMISES  ->  does not
#        predict the gate            (repetition vs ESMFold pLDDT, generated
#                                     sequences, r = +0.2125, n.s., n = 42)
#   (b)  the gate                                            ->  does not
#        predict the truth           (mean pLDDT vs real measured stability,
#                                     external designed set, rho = +0.205,
#                                     beaten by sequence length at +0.306)
#
# Splitting this into two figures destroys the entire reason it exists.
# CLAUDE.md section 8 names "reads as two half-papers" as a live risk.
#
# SOURCES
#   (a) analysis/figure1_scatter_data.csv  (42 rows; 8 of the 50 have an empty
#       generated field and are already excluded -- NEVER report this as n=50)
#       starred pair + correlation: paper-facts sections 5a / 5b
#   (b) paper-facts section 6 (Spearman table, n = 1,859)
#
# ############################################################################
# TWO LENGTHS PER SEQUENCE -- do not mix them (paper-facts section 5).
#   123 / 246 aa = the GENERATED CONTINUATION. The repetition metric is
#                  computed on this, so it is what panel (a) uses.
#   135 / 258 aa = the FOLDED CHAIN (prompt + continuation). Belongs to the
#                  structure panel in the appendix, NOT here.
# ############################################################################
#
# ############################################################################
# WHY PANEL (b) IS NOT A SCATTER.
#   Plotting measured stability against mean pLDDT and colouring by
#   design/scramble gives two cleanly separated clouds -- designs at 75.5
#   pLDDT and +0.289 stability, scrambles at ~63 and ~0. A reader concludes
#   pLDDT PREDICTS STABILITY, which is the opposite of the finding. The
#   separation is BETWEEN groups; there is nothing useful WITHIN them.
#   The ranked correlation is the honest form, and it needs no arithmetic.
# ############################################################################
#
# Output: fig3_gate.pdf. Run from the analysis/ folder.
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# palette -- matches Figure 1 (the draw.io pipeline diagram)
BLUE_FILL, BLUE_EDGE = "#DAE8FC", "#6C8EBF"   # standard pipeline
RED_FILL, RED_EDGE = "#F8CECC", "#B85450"     # the collapse verdict / the gate
GREY_BAR = "#CBCBCB"
GREY = "#8A8A8A"
DARK = "#222222"

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
})

fig = plt.figure(figsize=(6.6, 3.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 0.88], wspace=0.36)
ax_sc = fig.add_subplot(gs[0, 0])
ax_rho = fig.add_subplot(gs[0, 1])

# ============================================================== panel (a)
rows = list(csv.DictReader(open("figure1_scatter_data.csv")))
lens = [float(r["gen_len"]) for r in rows]
lmin, lmax = min(lens), max(lens)

# Marker area encodes generated length -- the mediator section 3 is about.
# This is a plotting scale, not a reported quantity.
#
# SQRT, not linear: gen_len spans 3 to 246 residues in this pool. A linear
# area map is dominated by the single 246-aa point and renders everything
# else as identical dots. Sqrt spreads the low end where most of the data is.
#
# NOTE the short tail is REAL and will be visible: several sequences are
# 3-18 aa, where distinct-3-mer is 1.0 more or less by construction (a 3-mer
# sequence contains exactly one 3-mer). That is not a plotting artifact --
# it is the metric being degenerate at short length, which is precisely the
# mediator argument section 3 makes. Do not hide it.
def area(n):
    return 12.0 + 75.0 * ((n - lmin) / (lmax - lmin)) ** 0.5

star = {}
for r in rows:
    x, y = float(r["distinct_3mer"]), float(r["plddt"])
    a = area(float(r["gen_len"]))
    collapsed = int(r["collapse"])
    ax_sc.scatter([x], [y], s=a,
                  facecolor=RED_FILL if collapsed else BLUE_FILL,
                  edgecolor=RED_EDGE if collapsed else BLUE_EDGE,
                  linewidth=0.8, alpha=0.92, zorder=3)
    if r["highlight"] in ("eerhv_repeat", "ubiquitin_repeat"):
        star[r["highlight"]] = (x, y, a)

ax_sc.axhline(60, ls=":", lw=0.9, color=GREY, zorder=1)
# label on the LEFT -- the right side is where the dense cluster sits
ax_sc.text(0.32, 58, "collapse threshold (pLDDT 60)", ha="left",
           va="top", fontsize=5.5, color=GREY)

# the matched pair: the two MOST repetitive sequences in the set are the
# worst and the best fold in it (paper-facts 5b)
for key in ("eerhv_repeat", "ubiquitin_repeat"):
    x, y, a = star[key]
    ax_sc.scatter([x], [y], s=a * 2.6, marker="*", facecolor="none",
                  edgecolor=DARK, linewidth=1.0, zorder=5)

ax_sc.annotate("0.339 $\\rightarrow$ 28.85\n123 aa", star["eerhv_repeat"][:2],
               textcoords="offset points", xytext=(12, -3), fontsize=5.6,
               color=DARK, linespacing=1.3, va="top")
ax_sc.annotate("0.365 $\\rightarrow$ 85.87\n246 aa", star["ubiquitin_repeat"][:2],
               textcoords="offset points", xytext=(12, 3), fontsize=5.6,
               color=DARK, linespacing=1.3, va="bottom")

# stats note, left-aligned in the empty bottom-left corner (nothing below y=20)
ax_sc.text(0.32, 3,
           "$r$ = +0.2125,  $t$(40) = +1.375,  n.s.  ($n$ = 42)\n"
           "marker area $\\propto$ generated length",
           ha="left", va="bottom", fontsize=5.4, color=DARK, linespacing=1.5)

ax_sc.set_xlim(0.30, 1.04)
ax_sc.set_ylim(0, 100)
ax_sc.set_xlabel("distinct 3-mer fraction  (less repetitive $\\rightarrow$)",
                 fontsize=7.2)
ax_sc.set_ylabel("ESMFold mean pLDDT", fontsize=7.2)
ax_sc.tick_params(labelsize=6.2)
ax_sc.set_title("what the steering optimises\ndoes not predict the gate",
                fontsize=7.6, color=DARK, pad=7, linespacing=1.4)

# ============================================================== panel (b)
# paper-facts section 6. Spearman with measured stability, n = 1,859.
feats = [
    # label,                rho,   is_the_gate
    ('"is it a design"',   0.313, False),
    ("sequence length",    0.306, False),
    ("pTM",                0.246, False),
    ("mean pLDDT",         0.205, True),
]

ys = list(range(len(feats)))[::-1]
for y, (label, rho, is_gate) in zip(ys, feats):
    ax_rho.barh(y, rho, height=0.62,
                color=RED_FILL if is_gate else GREY_BAR,
                edgecolor=RED_EDGE if is_gate else DARK,
                linewidth=0.7, zorder=2)
    ax_rho.text(rho + 0.008, y, "+%.3f" % rho, va="center", ha="left",
                fontsize=6.0, color=RED_EDGE if is_gate else DARK)

ax_rho.set_yticks(ys)
ax_rho.set_yticklabels([f[0] for f in feats], fontsize=6.4)
ax_rho.get_yticklabels()[-1].set_color(RED_EDGE)
ax_rho.set_xlim(0, 0.37)
ax_rho.set_ylim(-0.9, len(feats) - 0.3)
ax_rho.set_xlabel("Spearman $\\rho$ with measured stability", fontsize=7.2)
ax_rho.tick_params(axis="x", labelsize=6.2)
ax_rho.tick_params(axis="y", length=0)
for side in ("right", "top"):
    ax_rho.spines[side].set_visible(False)

ax_rho.text(0.365, -0.86, "$n$ = 1,859 external designed sequences",
            ha="right", va="bottom", fontsize=5.4, color=GREY)

ax_rho.set_title("and the gate does not\npredict the truth",
                 fontsize=7.6, color=DARK, pad=7, linespacing=1.4)

# panel letters
for ax, letter in ((ax_sc, "a"), (ax_rho, "b")):
    ax.text(-0.14, 1.17, "(%s)" % letter, transform=ax.transAxes,
            fontsize=8, fontweight="bold", color=DARK, va="top")

fig.savefig("fig3_gate.pdf", bbox_inches="tight")
print("wrote fig3_gate.pdf")
