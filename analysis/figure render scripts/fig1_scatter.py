# fig1_scatter.py -- Figure 1(a): "Repetition versus structural confidence"
#
# Reads analysis/figure1_scatter_data.csv (run from the analysis/ folder, or copy
# the CSV next to this script). Every value plotted comes from that CSV; the two
# emphasized points are the matched pair tagged in the `highlight` column.
# Output: fig1_scatter.pdf  (sized for a 0.48\linewidth panel).
import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("figure1_scatter_data.csv")))
hx, hy, cx, cy, pair = [], [], [], [], {}
for r in rows:
    x, y = float(r["distinct_3mer"]), float(r["plddt"])
    if r["highlight"] in ("eerhv_repeat", "ubiquitin_repeat"):
        pair[r["highlight"]] = (x, y)
    if int(r["collapse"]):
        cx.append(x); cy.append(y)
    else:
        hx.append(x); hy.append(y)

fig, ax = plt.subplots(figsize=(3.2, 3.0))
ax.axhline(60, ls=":", lw=1, color="0.5", zorder=1)
ax.text(1.01, 61, "collapse threshold (pLDDT 60)", ha="right", va="bottom",
        fontsize=6.5, color="0.4")
ax.scatter(hx, hy, s=26, facecolor="#2c7fb8", edgecolor="none", alpha=0.85,
           label="folds (pLDDT >= 60)", zorder=2)
ax.scatter(cx, cy, s=26, facecolor="none", edgecolor="#e34a33", linewidth=1.1,
           label="collapses (pLDDT < 60)", zorder=2)
for x, y in pair.values():
    ax.scatter([x], [y], s=150, marker="*", facecolor="gold",
               edgecolor="black", linewidth=0.8, zorder=5)
ax.annotate("0.339 -> 28.85", pair["eerhv_repeat"], textcoords="offset points",
            xytext=(9, -3), fontsize=7)
ax.annotate("0.365 -> 85.87", pair["ubiquitin_repeat"], textcoords="offset points",
            xytext=(9, 3), fontsize=7)
ax.set_xlabel("distinct 3-mer fraction (repetition)", fontsize=9)
ax.set_ylabel("ESMFold pLDDT (structural confidence)", fontsize=9)
ax.set_xlim(0.30, 1.03); ax.set_ylim(0, 100); ax.tick_params(labelsize=8)
ax.legend(fontsize=6.5, loc="lower right", framealpha=0.9)
fig.tight_layout(pad=0.4)
fig.savefig("fig1_scatter.pdf", bbox_inches="tight")
print("wrote fig1_scatter.pdf")
