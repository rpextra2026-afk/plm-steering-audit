# fig_appendix.py -- Appendix figures A2 to A5.
#
# A1 (the two structures at one shared zoom scale) is NOT here. It is a PyMOL
# render and its PDF already exists at analysis/final/fig1_structures.pdf.
# Do not delete that file -- it is the only rendered copy, and remaking it
# needs PyMOL plus analysis/fig1/fig1_{ubiquitin,eerhv}_repeat.pdb.
#
#   A2  sparse features vs raw activations at matched budget   (a REPLICATION)
#   A3  the design-vs-scramble gap sweep                       (UNDERPOWERED)
#   A4  the paired truncation                                  (the mediator)
#   A5  within-model layer dose asymmetry                      (why a claim died)
#
# SELF-CONTAINED. Every value is transcribed from paper/paper-facts.md:
#   A2  section 3b        A3  section 6 (gap sweep table)
#   A4  section 5c        A5  section 2e (asymmetry table)
#
# TWO NUMBERS DELIBERATELY LEFT FOR THE CAPTIONS, not drawn:
#   A2  dictionary health -- FVU 0.093, mean L0 32.0, 46.7% dead, 2,560
#       latents, 721k training tokens. Also the mandatory caveat: the gap
#       HALVED (-0.117 -> -0.049) as the dictionary improved, so this is an
#       UPPER BOUND on the disadvantage of sparse features, not a floor.
#   A3  at gap 1.25 the activations CI is [50.0%, 69.2%] -- lower bound
#       exactly at chance. Notebook verdict: real_trend_underpowered.
#
# Writes four separate PDFs so each can be placed independently.
# Run from anywhere -- no input files.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# palette -- matches Figure 1 and the body figures
TEAL, TEAL_FILL = "#2E7373", "#DCEDED"
RED_FILL, RED_EDGE = "#F8CECC", "#B85450"
BLUE_FILL, BLUE_EDGE = "#DAE8FC", "#6C8EBF"
PURPLE, PURPLE_FILL = "#3333AA", "#E4E4F4"
GREY_BAR, GREY, GREY_MID, DARK = "#CBCBCB", "#8A8A8A", "#9E9E9E", "#222222"

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
})


def finish(fig, ax_list, name):
    for ax in ax_list:
        for side in ("right", "top"):
            ax.spines[side].set_visible(False)
    fig.savefig(name, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print("wrote " + name)


# ==========================================================================
# A2 -- sparse features vs raw activations, matched budget, layer 30
# paper-facts 3b. THIS IS A REPLICATION of established prior work in language
# models. The direction is NOT ours. The caption must say so.
# Layer-12 rows are NOT plotted: that dictionary is 94.8% dead and its
# comparison is inadmissible (paper-facts 0.2 item 7).
# ==========================================================================
fig = plt.figure(figsize=(3.3, 2.6))
ax = fig.add_subplot(111)

ks = [3, 20, 64]
sae = [0.643, 0.749, 0.789]
raw = [0.683, 0.834, 0.851]
xs = [0, 1, 2]

ax.fill_between(xs, sae, raw, color=TEAL_FILL, alpha=0.6, lw=0, zorder=0)
ax.plot(xs, raw, color=TEAL, lw=1.8, marker="o", ms=4.2, label="raw activations")
ax.plot(xs, sae, color=GREY_MID, lw=1.5, marker="s", ms=3.8, ls="--",
        label="SAE latents")

for x, s_, r_ in zip(xs, sae, raw):
    ax.text(x, r_ + 0.011, "%.3f" % r_, ha="center", va="bottom",
            fontsize=5.4, color=TEAL)
    ax.text(x, s_ - 0.013, "%.3f" % s_, ha="center", va="top",
            fontsize=5.4, color=GREY_MID)

ax.set_xticks(xs)
ax.set_xticklabels([str(k) for k in ks], fontsize=6.4)
ax.set_xlim(-0.35, 2.5)
ax.set_ylim(0.60, 0.90)
ax.set_xlabel("feature budget $k$", fontsize=7.2)
ax.set_ylabel("collapse-prediction AUC", fontsize=7.2)
ax.tick_params(axis="y", labelsize=6.2)
# legend top-left: the only region with no data at any budget
ax.legend(fontsize=5.8, loc="upper left", frameon=False)
# one short line only -- dictionary health goes in the caption
ax.text(-0.30, 0.604,
        "layer 30  |  matched budget and classifier  |  in-fold selection",
        ha="left", va="bottom", fontsize=5.2, color=GREY)
ax.set_title("sparse features underperform raw,\nat every budget tested",
             fontsize=7.4, color=DARK, pad=6, linespacing=1.4)
finish(fig, [ax], "figA2_sae_vs_raw.pdf")


# ==========================================================================
# A3 -- the design-vs-scramble gap sweep
# paper-facts 6. MAY NOT be drawn as an established trend. The chance line is
# mandatory, and the collapsing pair count is the reason it is underpowered.
# ==========================================================================
fig = plt.figure(figsize=(3.5, 2.7))
ax = fig.add_subplot(111)

gaps = [0.25, 0.50, 0.75, 1.00, 1.25]
pairs = [724, 433, 253, 134, 71]
arm_a = [50.1, 48.0, 47.2, 46.0, 49.9]
arm_b = [51.5, 51.4, 52.4, 53.9, 58.5]
arm_c = [50.6, 49.9, 49.9, 51.2, 54.4]

ax.axhline(50, color=DARK, ls=":", lw=1.0, zorder=1)
# label parked to the right of the last data point -- it sat on the lines
ax.text(1.30, 50.2, "chance", ha="left", va="bottom", fontsize=5.5, color=DARK)

ax.plot(gaps, arm_b, color=TEAL, lw=1.7, marker="o", ms=3.8, label="activations")
ax.plot(gaps, arm_c, color=GREY_MID, lw=1.3, marker="D", ms=3.0, label="both")
ax.plot(gaps, arm_a, color=RED_EDGE, lw=1.3, marker="^", ms=3.4, ls="--",
        label="pLDDT features")

# the collapsing sample size is the whole reason this is underpowered
for g, n in zip(gaps, pairs):
    ax.text(g, 44.6, "n=%d" % n, ha="center", va="bottom",
            fontsize=5.0, color=GREY)

ax.set_xlim(0.17, 1.52)
ax.set_ylim(44, 66)
ax.set_xticks(gaps)
ax.set_xticklabels(["0.25", "0.50", "0.75", "1.00", "1.25"], fontsize=6.2)
ax.set_xlabel("minimum measured stability gap within a pair", fontsize=7.2)
ax.set_ylabel("balanced accuracy (%)", fontsize=7.2)
ax.tick_params(axis="y", labelsize=6.2)
ax.legend(fontsize=5.8, loc="upper left", frameon=False)
ax.set_title("coarse differences are seen, fine ones missed\n"
             "(underpowered: no level clears chance alone)",
             fontsize=7.2, color=DARK, pad=6, linespacing=1.4)
finish(fig, [ax], "figA3_gap_sweep.pdf")


# ==========================================================================
# A4 -- the paired truncation
# paper-facts 5c. THE SAME 50 SEQUENCES truncated to each length -- a paired
# within-sequence design. That pairing is what distinguishes this from the
# 2022 AlphaFold work, which varied length across DIFFERENT random sequences.
# The phenomenon is theirs; the paired design and the false-null consequence
# are ours.
# ==========================================================================
# SINGLE AXIS ONLY. The collapse-rate series was removed 2026-09-05: it is
# mean pLDDT thresholded at 60, i.e. the SAME measurement plotted twice, and
# on a twin axis whose range made the two series cross. That invited a
# reader to see a relationship the data does not assert. One line carries
# the claim; the annotation carries the statistics.
fig = plt.figure(figsize=(3.5, 2.7))
ax = fig.add_subplot(111)

trunc = [25, 40, 60, 80]
plddt = [61.85, 58.94, 55.91, 56.65]

ax.plot(trunc, plddt, color=BLUE_EDGE, lw=1.8, marker="o", ms=4.2, zorder=3)

for t, p in zip(trunc, plddt):
    ax.text(t, p + 0.40, "%.2f" % p, ha="center", va="bottom",
            fontsize=5.4, color=BLUE_EDGE)

ax.set_xlim(18, 87)
ax.set_ylim(51.5, 65)
ax.set_xticks(trunc)
ax.set_xticklabels([str(t) for t in trunc], fontsize=6.2)
ax.set_xlabel("truncation length (residues)", fontsize=7.2)
ax.set_ylabel("mean pLDDT", fontsize=7.2, color=BLUE_EDGE)
ax.tick_params(axis="y", labelsize=6.2, colors=BLUE_EDGE)

ax.text(19.5, 51.9,
        "the SAME 50 sequences at each length\n"
        "paired +5.21 pLDDT short-vs-long\n"
        "Wilcoxon $p$ = 0.0162  |  32/50 higher when truncated",
        ha="left", va="bottom", fontsize=5.2, color=DARK, linespacing=1.5)
ax.set_title("shorter sequences score higher\non the same metric",
             fontsize=7.4, color=DARK, pad=6, linespacing=1.4)
finish(fig, [ax], "figA4_truncation.pdf")


# ==========================================================================
# A5 -- within-model layer dose asymmetry
# paper-facts 2e. This is a MEASUREMENT ARTIFACT found in our own earlier
# analysis, and the reason the early-layers-are-more-sensitive claim died
# (blocklist item 2). Never present it as a finding about layers.
# ==========================================================================
fig = plt.figure(figsize=(3.7, 2.6))
ax = fig.add_subplot(111)

rows = [
    # model,          early alpha_rel, late alpha_rel, ratio, early/late layers
    ("Mistral-Prot",   0.010,  0.002, "4.47x", "2 / 7"),
    ("p-IgGen",       95.079, 26.774, "3.55x", "1 / 3"),
    ("RITA",          24.652, 13.572, "1.82x", "3 / 11"),
    ("ZymCTRL",       11.419,  7.065, "1.62x", "12 / 30"),
    ("ProtGPT2",       0.200,  0.156, "1.28x", "12 / 30"),
]

for i, (name, early, late, ratio, layers) in enumerate(rows):
    is_main = (name == "ProtGPT2")
    col = PURPLE if is_main else DARK
    ax.plot([late, early], [i, i], color=col, lw=1.1, zorder=2)
    # ProtGPT2's ratio is only 1.28x, so on a log axis its two markers very
    # nearly coincide and it reads as having a single point while every other
    # model shows two. Shrink both markers and drop the fill alpha slightly so
    # the pair stays distinguishable at print size.
    tight = (early / late) < 1.5
    s_late, s_early = (17, 20) if tight else (26, 30)
    ax.scatter([late], [i], s=s_late, facecolor="white", edgecolor=col,
               linewidth=1.1, zorder=4)
    ax.scatter([early], [i], s=s_early, facecolor=col, edgecolor=col,
               linewidth=1.1, zorder=3)
    ax.text(1.7e2, i, ratio, ha="left", va="center", fontsize=5.8, color=col)

ax.set_yticks(range(len(rows)))
ax.set_yticklabels(["%s\n(%s)" % (r[0], r[4]) for r in rows], fontsize=5.8,
                   linespacing=1.3)
ax.get_yticklabels()[-1].set_color(PURPLE)
ax.set_xscale("log")
ax.set_xlim(8e-4, 6.0e2)
ax.set_ylim(-0.75, len(rows) - 0.25)
ax.set_xlabel(r"relative push $\alpha_{rel}$ at nominal 1$\times$", fontsize=7.2)
ax.tick_params(axis="x", labelsize=6.0)
ax.tick_params(axis="y", length=0)

# legend goes UPPER LEFT -- it collided with the 4.47x label at lower right,
# and no model has data below alpha_rel 0.05 on the top row
ax.scatter([], [], s=30, facecolor=DARK, edgecolor=DARK, label="early layer")
ax.scatter([], [], s=26, facecolor="white", edgecolor=DARK, label="late layer")
ax.legend(fontsize=5.6, loc="upper left", frameon=False, scatterpoints=1,
          handletextpad=0.4, borderpad=0.2)

ax.set_title("every model receives more relative push early\n"
             "(median 1.82x), so matched push scales per layer",
             fontsize=7.0, color=DARK, pad=6, linespacing=1.4)
finish(fig, [ax], "figA5_layer_asymmetry.pdf")

print("")
print("A1 is a PyMOL render and is NOT produced here.")
print("It already exists at analysis/final/fig1_structures.pdf -- do not delete it.")
