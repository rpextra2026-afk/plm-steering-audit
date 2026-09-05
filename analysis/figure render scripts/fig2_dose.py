# fig2_dose.py -- Figure 2: dose-response cliff + prior-method operating points
#
# Self-contained: all values are transcribed from paper/paper-facts.md
#   - 8 dose points + Wilson CIs: paper-facts.md section 2a
#   - fitted floored logistic (floor 0.640, k 9.48, alpha50 1.240, CI [1.138, 1.321]): 2a
#   - prior-method bands (alpha=1 -> 0.563-0.753x, alpha=2 -> 1.126-1.506x): 2d
# Runs anywhere with matplotlib. Output: fig2_dose.pdf (sized for 0.6\linewidth).
import matplotlib, math
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dose  = [0.00, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00]
rate  = [64.0, 50.0, 50.0, 66.0, 82.0, 100.0, 98.0, 100.0]
ci_lo = [50.1, 36.6, 36.6, 52.2, 69.2, 92.9, 89.5, 92.9]
ci_hi = [75.9, 63.4, 63.4, 77.6, 90.2, 100.0, 99.6, 100.0]
yerr  = [[r - lo for r, lo in zip(rate, ci_lo)], [hi - r for r, hi in zip(rate, ci_hi)]]

floor, k, a50, a50ci = 0.640, 9.48, 1.240, (1.138, 1.321)
curve = lambda x: 100.0 * (floor + (1 - floor) / (1 + math.exp(-k * (x - a50))))
xs = [i * 2.15 / 200 for i in range(201)]
ys = [curve(x) for x in xs]

fig, ax = plt.subplots(figsize=(4.2, 2.8))
ax.axvspan(0.563, 0.753, color="#2ca25f", alpha=0.15, lw=0, label=r"prior $\alpha=1$")
ax.axvspan(1.126, 1.506, color="#f16913", alpha=0.15, lw=0, label=r"prior $\alpha=2$")
ax.axvspan(a50ci[0], a50ci[1], color="0.6", alpha=0.25, lw=0)
ax.axvline(a50, color="black", ls="--", lw=1, label=r"$\alpha_{50}=1.24\times$")
ax.axhline(64.0, color="0.5", ls=":", lw=1)
ax.plot(xs, ys, color="black", lw=1.8, label="floored logistic fit")
ax.errorbar(dose, rate, yerr=yerr, fmt="o", color="black", ms=4, capsize=2, lw=1,
            label="empirical (N=50/dose)")
ax.set_xlabel("steering dose (multiples of reference norm)", fontsize=9)
ax.set_ylabel("collapse rate (%)", fontsize=9)
ax.set_xlim(-0.03, 2.15); ax.set_ylim(0, 105); ax.tick_params(labelsize=8)
ax.legend(fontsize=6.5, loc="center right", framealpha=0.9)
fig.tight_layout(pad=0.4)
fig.savefig("fig2_dose.pdf", bbox_inches="tight")
print("wrote fig2_dose.pdf")
