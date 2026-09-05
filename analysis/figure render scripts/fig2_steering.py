# fig2_steering.py -- Figure 2: the steering audit, three panels.
#
# Reads top-to-bottom as one argument:
#   strip -> nominal dose is meaningless
#   (a)   -> here is the measured curve, and where the published method falls on it
#   (b)   -> direction only matters below the cliff
#
# SELF-CONTAINED. Every value below is transcribed from paper/paper-facts.md:
#   strip     section 2e  (NB34)  per-model alpha_rel, 41,033x spread
#   panel (a) section 2a  (NB37)  8 dose triples, Wilson CIs, floored-logistic fit
#             section 2b         the published method alpha=1 / alpha=2 bands
#   panel (b) section 2c  (NB35, ORTHOGONAL from the NB39 refold)  8 conditions
#
# ############################################################################
# THREE DIFFERENT CONTROLS APPEAR IN THIS FIGURE. THEY ARE NOT INTERCHANGEABLE.
#   64.0%  NB37 own control -> panel (a), dose 0.00x
#   52.0%  NB35 own control -> panel (b), and EVERY significance claim in it
#   62.0%  NB03 control     -> belongs to NEITHER panel. Never use it here.
# ############################################################################
#
# PALETTE -- matches Figure 1 (the draw.io pipeline diagram):
#   purple  = the published steering method   (real vector; ProtGPT2; alpha=1 band)
#   red     = the collapse verdict            (alpha=2 band -- it straddles the cliff)
#   grey    = everything else
#
# Output: fig2_steering.pdf  (full linewidth). Supersedes fig2_dose.pdf.
# Run from the analysis/ folder, then move the PDF into analysis/final/.
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PURPLE = "#3333AA"       # published steering method  (Fig 1: h + alpha*v, layer 12)
PURPLE_FILL = "#E4E4F4"  # Fig 1: layer-12 band tint
RED_FILL = "#F8CECC"     # Fig 1: collapse box
RED_EDGE = "#B85450"
GREY_BAR = "#CBCBCB"
GREY = "#8A8A8A"
DARK = "#222222"

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
})

fig = plt.figure(figsize=(6.6, 4.3))
gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 2.6],
                      width_ratios=[1.12, 0.88], hspace=0.80, wspace=0.30)
ax_strip = fig.add_subplot(gs[0, :])
ax_dose = fig.add_subplot(gs[1, 0])
ax_dir = fig.add_subplot(gs[1, 1])

# ---------------------------------------------------------------- strip
# Five segments, one per model, each spanning that model's two tested layers.
# Ranges are disjoint, so one row suffices. Labels ALTERNATE above/below --
# ZymCTRL, RITA and p-IgGen sit within ~0.4 decades and collide on one side.
models = [
    # name,           alpha_rel low, alpha_rel high, label above the line?
    ("Mistral-Prot",   0.002,  0.010, True),
    ("ProtGPT2",       0.156,  0.200, False),
    ("ZymCTRL",        7.065, 11.419, True),
    ("RITA",          13.572, 24.652, False),
    ("p-IgGen",       26.774, 95.079, True),
]

for name, lo, hi, above in models:
    is_main = (name == "ProtGPT2")
    colour = PURPLE if is_main else DARK
    width = 2.8 if is_main else 1.4
    ax_strip.plot([lo, hi], [0, 0], color=colour, lw=width,
                  solid_capstyle="butt", zorder=3)
    for x in (lo, hi):
        ax_strip.plot([x, x], [-0.14, 0.14], color=colour,
                      lw=width * 0.65, zorder=3)
    ax_strip.text(math.sqrt(lo * hi), 0.30 if above else -0.40, name,
                  ha="center", va="bottom" if above else "top",
                  fontsize=6.2, color=colour,
                  fontweight="bold" if is_main else "normal")

# alpha_rel = 1.0 is a DEFINITION, not a result: the added vector equals the
# hidden state it is added to. The ProtGPT2 cliff is deliberately NOT drawn
# across this strip -- it is model-specific (paper-facts 2e).
ax_strip.axvline(1.0, color=GREY, ls=":", lw=0.9, zorder=1)
ax_strip.text(1.0, 0.62, r"$\|v\|=\|h\|$", ha="center", va="bottom",
              fontsize=5.8, color=GREY,
              bbox=dict(facecolor="white", edgecolor="none", pad=1.2))

# span bracket, with its label BELOW the bracket rather than sitting on it
b_lo, b_hi, b_y = 0.002, 95.079, -1.15
ax_strip.plot([b_lo, b_hi], [b_y, b_y], color=GREY, lw=0.8)
for x in (b_lo, b_hi):
    ax_strip.plot([x, x], [b_y, b_y + 0.16], color=GREY, lw=0.8)
ax_strip.text(math.sqrt(b_lo * b_hi), b_y - 0.16,
              r"41,033$\times$ spread, every condition labelled 1$\times$",
              ha="center", va="top", fontsize=6.2, color=DARK)

ax_strip.set_xscale("log")
ax_strip.set_xlim(0.0011, 300)
ax_strip.set_ylim(-2.0, 1.05)
ax_strip.set_yticks([])
for side in ("left", "right", "top"):
    ax_strip.spines[side].set_visible(False)
ax_strip.tick_params(axis="x", labelsize=6, pad=1)
ax_strip.set_xlabel(r"relative push $\alpha_{rel}$  (fraction of $\|h\|$ displaced)",
                    fontsize=6.8, labelpad=1)

# ------------------------------------------------------------- panel (a)
dose = [0.00, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00]
rate = [64.0, 50.0, 50.0, 66.0, 82.0, 100.0, 98.0, 100.0]
ci_lo = [50.1, 36.6, 36.6, 52.2, 69.2, 92.9, 89.5, 92.9]
ci_hi = [75.9, 63.4, 63.4, 77.6, 90.2, 100.0, 99.6, 100.0]
yerr = [[r - lo for r, lo in zip(rate, ci_lo)],
        [hi - r for r, hi in zip(rate, ci_hi)]]

FLOOR, K, A50 = 0.640, 9.48, 1.240
A50_CI = (1.138, 1.321)
curve = lambda x: 100.0 * (FLOOR + (1 - FLOOR) / (1 + math.exp(-K * (x - A50))))
xs = [i * 2.15 / 300 for i in range(301)]

# the published method operating points, in the Figure 1 palette:
# alpha=1 purple (the method itself), alpha=2 red (it straddles the cliff)
ax_dose.axvspan(0.563, 0.753, color=PURPLE_FILL, lw=0, zorder=0)
ax_dose.axvspan(1.126, 1.506, color=RED_FILL, lw=0, zorder=0)
ax_dose.text(0.658, 70, r"their $\alpha$=1", rotation=90, ha="center",
             va="bottom", fontsize=5.8, color=PURPLE)
ax_dose.text(1.316, 30, r"their $\alpha$=2", rotation=90, ha="center",
             va="bottom", fontsize=5.8, color=RED_EDGE)

ax_dose.axvspan(A50_CI[0], A50_CI[1], color="0.55", alpha=0.30, lw=0, zorder=1)
ax_dose.axvline(A50, color=DARK, ls="--", lw=1.1, zorder=2)
ax_dose.text(A50, 108, r"$\alpha_{50}=1.24\times$", ha="center", va="bottom",
             fontsize=6.2, color=DARK)

ax_dose.axhline(64.0, color=GREY, ls=":", lw=0.9, zorder=1)
ax_dose.text(2.12, 66, "NB37 control 64.0%", ha="right", va="bottom",
             fontsize=5.5, color=GREY)

ax_dose.plot(xs, [curve(x) for x in xs], color=DARK, lw=1.6, zorder=3)
ax_dose.errorbar(dose, rate, yerr=yerr, fmt="o", color=DARK, ms=3.4,
                 capsize=2.0, lw=1.0, zorder=4)

# the two significance annotations, as a compact note -- no leader lines,
# so nothing crosses the fitted curve or the shaded bands
ax_dose.text(0.03, 16,
             "1.25$\\times$  n.s. vs control ($p$ = 0.0705)\n"
             "1.50$\\times$  first significant ($p$ < 0.0001)",
             ha="left", va="bottom", fontsize=5.6, color=DARK, linespacing=1.5)

ax_dose.set_xlim(-0.04, 2.15)
ax_dose.set_ylim(0, 120)
ax_dose.set_yticks([0, 20, 40, 60, 80, 100])
ax_dose.set_xlabel(r"steering dose ($\times$ reference norm)", fontsize=7.2)
ax_dose.set_ylabel("collapse rate (%)", fontsize=7.2)
ax_dose.tick_params(labelsize=6.2)

# secondary axis -- ticks ONLY at locked (dose, alpha_rel) pairs from the
# section 2a table, so nothing is interpolated
sec = ax_dose.secondary_xaxis("top")
sec.set_xticks([0.00, 1.00, 2.00])
sec.set_xticklabels(["0.000", "0.172", "0.344"], fontsize=5.6)
sec.set_xlabel(r"$\alpha_{rel}$", fontsize=6.2, labelpad=2)
sec.tick_params(length=2, pad=1)

# ------------------------------------------------------------- panel (b)
# NB35, N=50/condition. CONTROL here is 52.0%, NB35 own -- NOT section 1 62.0%.
conds = [
    # x, label,        rate,  ci_lo, ci_hi, is_real_vector
    (0, "control",      52.0, 38.5, 65.2, False),
    (1, "real",         64.0, 50.1, 75.9, True),
    (2, "random",       92.0, 81.2, 96.8, False),
    (4, "real",        100.0, 92.9, 100.0, True),
    (5, "random",      100.0, 92.9, 100.0, False),
    (6, "shuffled",    100.0, 92.9, 100.0, False),
    (7, "orthogonal",  100.0, 92.9, 100.0, False),
    (8, "negated",      70.0, 56.2, 80.9, False),
]

for x, label, r, lo, hi, is_real in conds:
    ax_dir.bar(x, r, width=0.76, color=PURPLE if is_real else GREY_BAR,
               edgecolor=DARK, linewidth=0.6, zorder=2)
    ax_dir.errorbar(x, r, yerr=[[r - lo], [hi - r]], fmt="none",
                    ecolor=DARK, capsize=1.8, lw=0.9, zorder=3)

ax_dir.axhline(52.0, color=GREY, ls=":", lw=0.9, zorder=1)
ax_dir.text(8.55, 53.5, "control", ha="left", va="bottom",
            fontsize=5.5, color=GREY)

# divider + group labels INSIDE the axes, so they never touch the tick labels
ax_dir.axvline(3.0, color="#DDDDDD", lw=0.8, zorder=0)
ax_dir.text(1.0, 120, r"1$\times$", ha="center", fontsize=7.5, color=DARK)
ax_dir.text(6.0, 120, r"2$\times$", ha="center", fontsize=7.5, color=DARK)

# bracket 1: REAL vs RANDOM at 1x -- non-overlapping CIs
ax_dir.plot([1, 1, 2, 2], [99, 103, 103, 99], color=DARK, lw=0.7)
ax_dir.text(1.5, 104, "CIs do not\noverlap", ha="center", va="bottom",
            fontsize=5.2, color=DARK, linespacing=1.1)

# bracket 2: the four saturating conditions at 2x
ax_dir.plot([4, 4, 7, 7], [104, 108, 108, 104], color=DARK, lw=0.7)
ax_dir.text(5.5, 109, "Fisher $p$ = 1.0000", ha="center", va="bottom",
            fontsize=5.8, color=DARK)

ax_dir.set_xticks([c[0] for c in conds])
ax_dir.set_xticklabels([c[1] for c in conds], fontsize=6, rotation=90)
ax_dir.set_xlim(-0.7, 10.3)
ax_dir.set_ylim(0, 132)
ax_dir.set_yticks([0, 20, 40, 60, 80, 100])
ax_dir.set_ylabel("collapse rate (%)", fontsize=7.2)
ax_dir.tick_params(axis="y", labelsize=6.2)

# pad_inches added 2026-09-05: the rotated "orthogonal" tick label sits at the
# very bottom of panel (b) and was clipped at the figure edge. That label is
# one of the four saturating conditions carrying the "direction is irrelevant
# at saturation" claim, so it must be fully legible.
fig.savefig("fig2_steering.pdf", bbox_inches="tight", pad_inches=0.08)
print("wrote fig2_steering.pdf")
