# fig4_predictor.py -- Figure 4: the predictor, and its limit.
#
#   (a)  restricting to the calls the classifier is most confident about buys
#        real precision. The gap between the curve and the pool base rate is
#        the value of abstention.        (NB59, N=1200, Wilson 95% CIs)
#   (b)  but the signal is absent in one of three models. RITA-small sits at
#        chance at every prefix length.  (NB42 / NB44)
#
# The two panels are deliberately a positive and its limit, in one frame --
# the same shape as Figure 3. A figure showing only (a) would overstate the
# one genuinely positive result in the paper.
#
# SELF-CONTAINED. Every value is transcribed from paper/paper-facts.md:
#   panel (a) section 3c  (NB59)  6 coverage rows, Wilson CIs, base rate 52.2%
#   panel (b) section 3d  (NB42/NB44)  three models x four prefix lengths
#
# ############################################################################
# TWO DIFFERENT 30s -- paper-facts section 3d.
#   the probe reads LAYER 30 (of the 36 layers of ProtGPT2)
#   early abort retains 94.9% of AUC at the first 30 RESIDUES
# Panel (b) is about RESIDUES. Never label it with a bare "30".
# ############################################################################
#
# ############################################################################
# NEVER QUOTE RITA RETENTION. Its 97.7% is retention of a CHANCE-LEVEL AUC and
# is meaningless -- the source CSV records signal_real as no_at_chance. The
# flat line is the honest content. Retention percentages are only
# interpretable when full-length AUC is above chance (ProtGPT2, ZymCTRL).
# ############################################################################
#
# ############################################################################
# BANNED in this figure and its caption: calibrated, conformal,
# distribution-free (paper-facts 0.1). This is an ordinary random forest
# evaluated on out-of-fold probabilities. PERMITTED and encouraged:
# selective prediction, abstention, risk-coverage, triage, early abort.
# ############################################################################
#
# Output: fig4_predictor.pdf. Run from anywhere -- no input files.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# palette -- matches Figure 1 (the draw.io pipeline diagram)
TEAL = "#2E7373"        # added by this study -- the triage classifier
TEAL_FILL = "#DCEDED"
RED_EDGE = "#B85450"    # the null / the limit
GREY = "#8A8A8A"
GREY_MID = "#9E9E9E"
DARK = "#222222"

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
})

fig = plt.figure(figsize=(6.6, 3.0))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.30)
ax_rc = fig.add_subplot(gs[0, 0])
ax_ea = fig.add_subplot(gs[0, 1])

# ============================================================== panel (a)
# paper-facts section 3c. NB59, N=1200, 10 repeats of stratified 10-fold.
cov = [10, 20, 30, 50, 70, 100]
prec = [98.1, 95.0, 92.6, 85.8, 79.4, 72.6]
lo = [87.9, 89.3, 87.7, 81.5, 75.5, 69.2]
hi = [99.6, 97.6, 95.3, 89.0, 82.7, 75.9]
BASE = 52.2
yerr = [[p - a for p, a in zip(prec, lo)], [b - p for p, b in zip(prec, hi)]]

# the gap between the curve and the base rate IS the value of abstention
ax_rc.fill_between(cov, [BASE] * len(cov), prec, color=TEAL_FILL,
                   alpha=0.75, lw=0, zorder=1)
ax_rc.axhline(BASE, color=GREY, ls=":", lw=0.9, zorder=2)
ax_rc.text(102, BASE + 1.2, "pool base rate 52.2%", ha="right", va="bottom",
           fontsize=5.5, color=GREY)

ax_rc.plot(cov, prec, color=TEAL, lw=1.6, zorder=3)
ax_rc.errorbar(cov, prec, yerr=yerr, fmt="o", color=TEAL, ms=3.6,
               capsize=2.0, lw=1.0, zorder=4)

ax_rc.text(103, 103, "top decile:  98.1%  [87.9, 99.6]  on 43 calls",
           ha="right", va="top", fontsize=5.6, color=DARK)

ax_rc.set_xlim(0, 105)
ax_rc.set_ylim(45, 104)
ax_rc.set_xticks([10, 20, 30, 50, 70, 100])
ax_rc.set_xlabel("coverage (% of the pool screened)", fontsize=7.2)
ax_rc.set_ylabel("precision on predicted-collapse calls (%)", fontsize=7.2)
ax_rc.tick_params(labelsize=6.2)
for side in ("right", "top"):
    ax_rc.spines[side].set_visible(False)
ax_rc.set_title("abstaining on the uncertain calls\nbuys precision",
                fontsize=7.6, color=DARK, pad=7, linespacing=1.4)

# ============================================================== panel (b)
# paper-facts section 3d. Three models, four prefix lengths.
# x is CATEGORICAL -- "full" is not a residue count.
xs = [0, 1, 2, 3]
xlabels = ["first 10", "first 20", "first 30", "full"]

protgpt2 = [0.682, 0.720, 0.740, 0.780]
pg_lo = [0.606, 0.637, 0.661, 0.709]
pg_hi = [0.734, 0.767, 0.784, 0.835]
zymctrl = [0.597, 0.643, 0.646, 0.658]
rita = [0.503, 0.496, 0.495, 0.506]

ax_ea.axhline(0.5, color=GREY, ls=":", lw=0.9, zorder=1)
ax_ea.text(-0.15, 0.492, "chance", ha="left", va="top",
           fontsize=5.5, color=GREY)

# min-max band is available for ProtGPT2 only (15 repeats per truncation)
ax_ea.fill_between(xs, pg_lo, pg_hi, color=TEAL_FILL, alpha=0.8, lw=0, zorder=2)
ax_ea.plot(xs, protgpt2, color=TEAL, lw=1.8, marker="o", ms=3.6, zorder=5)
ax_ea.plot(xs, zymctrl, color=GREY_MID, lw=1.4, marker="s", ms=3.2, zorder=4)
ax_ea.plot(xs, rita, color=RED_EDGE, lw=1.4, marker="^", ms=3.4,
           ls="--", zorder=4)

# RITA sits on the chance line, so its label drops below the marker
# rather than being centred on it
for y, name, colour, va in ((protgpt2[-1], "ProtGPT2", TEAL, "center"),
                           (zymctrl[-1], "ZymCTRL", GREY_MID, "center"),
                           (rita[-1] - 0.028, "RITA-small", RED_EDGE, "top")):
    ax_ea.text(3.12, y, name, ha="left", va=va,
               fontsize=6.2, color=colour)

ax_ea.set_xlim(-0.20, 4.45)
ax_ea.set_ylim(0.44, 0.87)
ax_ea.set_xticks(xs)
ax_ea.set_xticklabels(xlabels, fontsize=6.2)
ax_ea.set_xlabel("residues of the sequence seen", fontsize=7.2)
ax_ea.set_ylabel("collapse-prediction AUC", fontsize=7.2)
ax_ea.tick_params(axis="y", labelsize=6.2)
for side in ("right", "top"):
    ax_ea.spines[side].set_visible(False)
ax_ea.set_title("but the signal is absent\nin one of three models",
                fontsize=7.6, color=DARK, pad=7, linespacing=1.4)

# panel letters
for ax, letter in ((ax_rc, "a"), (ax_ea, "b")):
    ax.text(-0.15, 1.17, "(%s)" % letter, transform=ax.transAxes,
            fontsize=8, fontweight="bold", color=DARK, va="top")

fig.savefig("fig4_predictor.pdf", bbox_inches="tight")
print("wrote fig4_predictor.pdf")
