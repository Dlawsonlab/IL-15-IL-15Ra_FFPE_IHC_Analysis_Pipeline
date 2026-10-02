#!/usr/bin/env python3
r"""
Build the final cohort figure. It plots IL-15 (panel F) and IL-15Ra (panel G) positive area on
a LINEAR % axis, uses the FOV-mean endpoint, and overlays every section on the bars. It
follows the geometry of panels F and G in 8_Biorxiv_submission/1_Figures/Fig_3_V6.ai, so it
can be placed in Illustrator at 100 %.

Two versions, identical except for labelling:
  fig_cohort_bars_pct_legend    x axis "Brain" / "Liver", state legend box at the right --
                                the reference layout exactly
  fig_cohort_bars_pct_nolegend  no legend; each bar carries its group name (rotated 45 deg,
                                because four names do not fit a 90 pt axis horizontally).
                                Plot area and bar positions are identical to the legend
                                version, so either drops into the same slot.

GEOMETRY is measured from Fig_3_V6.ai (pdftotext -bbox, and a 600 dpi render):
  plot area 90 x 90 pt; bar width 10.4 pt; bar centres 11.2 / 27.9 / 61.5 / 78.3 pt from the
  y axis (16.7 pt within an organ pair); all strokes 0.75 pt; dots 4.3 pt across; ticks 3.7 pt;
  legend box bottom on the x axis, 4.1 pt right of the plot; Arial 7 pt (ticks, labels,
  legend), 7 pt bold (panel titles), 12 pt bold (panel letters). Brackets stack from 3.1 pt
  above the plot in ~6 pt steps. Only the lowest level carries end drops, as in the reference.
  The p-values are 5 pt, a deliberate departure from the reference. The reference's 4.3 pt is
  below the minimum type size of most journals.

Y AXIS starts at 0 and is capped at the highest section plus one dot radius. That dot then
sits inside the axis and does not hang half over its end, and the bars use the full height.
There is no tick at the cap, so ticks stay on round values.

STATISTICS follow METHODS section 10.
  Kruskal-Wallis runs across the four groups, and a two-sided Mann-Whitney U runs for all six
  group pairs. Benjamini-Hochberg FDR corrects across those six. The correction family is the
  six tests PERFORMED. Contrasts reaching adjusted p < 0.05 are drawn, and the rest are
  omitted but still corrected for. Bars are group MEDIANS and whiskers the INTERQUARTILE
  RANGE, because the tests are rank-based (METHODS section 10 states this). The linear axis
  changes only how the data look. Every p-value here equals those in cohort_stats.txt /
  cohort_stats_fov_dots.txt.

Dots are placed by a deterministic constrained beeswarm with no randomness. The figure is
therefore byte-identical run to run, and no two dots overlap unless the bar width forces it.

TYPE stays editable in Illustrator. pdf/ps fonttype 42 and svg fonttype 'none' keep every
label as live Arial text. Mathtext is pinned to Arial so the superscript in "IL-15+" does not
fall back to DejaVu Sans.

Usage:  python3 scripts/figures/plot_cohort_bars_pct.py
Writes: outputs/figures/cohort/fig_cohort_bars_pct_legend.{pdf,svg,png}
        outputs/figures/cohort/fig_cohort_bars_pct_nolegend.{pdf,svg,png}
        outputs/tables/cohort_stats_bars_pct.txt     every number the figure prints
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(ROOT)

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["mathtext.fontset"] = "custom"
for _k in ("rm", "it", "bf", "sf"):
    matplotlib.rcParams[f"mathtext.{_k}"] = "Arial"
matplotlib.rcParams["mathtext.default"] = "regular"

# --- reference palette and geometry (points) ------------------------------------------------
RED, RED_FILL, GREY_FILL, INK = "#ED2024", "#F57E80", "#7F7F7F", "#000000"
LW = 0.75
AX_W = AX_H = 90.0
XC = [11.2, 27.9, 61.5, 78.3]           # bar centres, pt from the y axis
BAR_W = 10.4
DIVIDER = (XC[1] + XC[2]) / 2
DOT_D = 4.3
TICK_LEN = 3.7
CAP_W = BAR_W / 2
LEVEL0, STEP, DROP = 3.1, 6.3, 2.6      # bracket stack above the plot
P_PT, TXT_PT, TITLE_PT, LETTER_PT = 5.0, 7.0, 7.0, 12.0

# (csv tissue_type, full group name, is-metastatic)
CELLS = [("Normal brain", "Non-neoplastic brain", False),
         ("BCBM",         "BCBM",                 True),
         ("Normal liver", "Non-neoplastic liver", False),
         ("Liver met",    "Liver metastasis",     True)]
PANELS = [("F", "fov_mean_brown_af",   r"IL-15$^{+}\;$area",
           "% IL-15 expression\n[IL-15 area/total area]"),
          ("G", "fov_mean_magenta_af", r"IL-15Rα$^{+}\;$area",
           "% IL-15Rα expression\n[IL-15Rα area/total area]")]
SUP = r"$^{+}\;$"
ALL_PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

rows = list(csv.DictReader(open("outputs/tables/cohort_table.csv")))


def fmt_p(p):
    return "<0.0001" if p < 0.0001 else f"{p:.4f}"


def bh(pvals):
    m = len(pvals); order = np.argsort(pvals); adj = [0.0] * m; prev = 1.0
    for k, idx in enumerate(reversed(order)):
        prev = min(prev, pvals[idx] * m / (m - k)); adj[idx] = prev
    return adj


def swarm(y_pt, half_w=BAR_W / 2, d=DOT_D):
    """Place dots by a constrained beeswarm in point units. Dots go in lowest value first. Each
    takes the smallest offset within +/- half_w that clears every dot already placed. If none
    does, it takes the offset with the most clearance. The result is deterministic, so the
    figure is reproducible."""
    cand = [0.0]
    k = 1
    while k * d / 2 <= half_w + 1e-9:
        cand += [k * d / 2, -k * d / 2]
        k += 1
    xs = np.zeros(len(y_pt)); placed = []
    for i in np.argsort(y_pt, kind="stable"):
        best, best_gap = cand[0], -1.0
        for c in cand:
            gap = min((np.hypot(c - xs[j], y_pt[i] - y_pt[j]) for j in placed), default=np.inf)
            if gap >= 0.92 * d:
                best = c
                break
            if gap > best_gap:
                best, best_gap = c, gap
        xs[i] = best
        placed.append(i)
    return xs


def analyse(col):
    vals = [np.array([float(r[col]) * 100 for r in rows if r["tissue_type"] == tt])
            for tt, _, _ in CELLS]
    H, p_kw = stats.kruskal(*vals)
    raw = [stats.mannwhitneyu(vals[i], vals[j], alternative="two-sided").pvalue
           for i, j in ALL_PAIRS]
    return vals, H, p_kw, raw, dict(zip(ALL_PAIRS, bh(raw)))


def stack(padj):
    """Pack significant brackets into levels. The shortest span goes lowest. Two brackets share
    a level when they do not overlap, so the within-organ pairs sit side by side."""
    levels, out = [], []
    for i, j in sorted((q for q in ALL_PAIRS if padj[q] < 0.05),
                       key=lambda q: (XC[q[1]] - XC[q[0]], q[0])):
        x0, x1 = XC[i], XC[j]
        for L, occ in enumerate(levels):
            if all(x1 + 3 < a or x0 - 3 > b for a, b in occ):
                occ.append((x0, x1)); out.append(((i, j), L)); break
        else:
            levels.append([(x0, x1)]); out.append(((i, j), len(levels) - 1))
    return out


RESULTS = {col: analyse(col) for _, col, _, _ in PANELS}
N_LEVELS = max(max((L for _, L in stack(RESULTS[c][4])), default=-1) + 1 for _, c, _, _ in PANELS)
TITLE_D = LEVEL0 + (N_LEVELS - 1) * STEP + 0.6 + 3.6 + 4.5   # title baseline above plot, pt


def build(with_legend):
    left, right = 62.0, (80.0 if with_legend else 20.0)
    above = TITLE_D + 10.0
    below = 16.0 if with_legend else 62.0
    gap = 14.0
    fw = left + AX_W + right
    fh = 6 + 2 * (above + AX_H + below) + gap + 6
    fig = plt.figure(figsize=(fw / 72, fh / 72))

    for pi, (letter, col, title, ylab) in enumerate(PANELS):
        vals, H, p_kw, raw, padj = RESULTS[col]
        bottom = fh - 6 - above - AX_H - pi * (above + AX_H + below + gap)
        ax = fig.add_axes([left / fw, bottom / fh, AX_W / fw, AX_H / fh])

        vmax = max(v.max() for v in vals)
        cap = vmax / (1 - (DOT_D / 2) / AX_H)       # top dot sits just inside the axis
        ax.set_xlim(0, AX_W)
        ax.set_ylim(0, cap)

        for x, v, (_, _, met) in zip(XC, vals, CELLS):
            edge = RED if met else INK
            med = np.median(v); q1, q3 = np.percentile(v, [25, 75])
            ax.bar(x, med, width=BAR_W, bottom=0, facecolor=RED_FILL if met else GREY_FILL,
                   edgecolor=edge, linewidth=LW, zorder=2)
            off = swarm(v / cap * AX_H)
            ax.scatter(x + off, v, s=DOT_D ** 2, facecolor=edge, edgecolor="none",
                       zorder=3, clip_on=False)
            ax.plot([x, x], [q1, q3], color=edge, lw=LW, zorder=4, solid_capstyle="butt")
            for q in (q1, q3):
                ax.plot([x - CAP_W / 2, x + CAP_W / 2], [q, q], color=edge, lw=LW, zorder=4,
                        solid_capstyle="butt")

        ax.axvline(DIVIDER, color=INK, lw=LW, ls=(0, (1.4, 1.7)), zorder=1)

        # significance brackets, above the plot area
        tr = ax.get_xaxis_transform()
        for (i, j), L in stack(padj):
            y = 1 + (LEVEL0 + L * STEP) / AX_H
            ax.plot([XC[i], XC[j]], [y, y], transform=tr, color=INK, lw=LW,
                    clip_on=False, solid_capstyle="butt")
            if L == 0:
                for e in (XC[i], XC[j]):
                    ax.plot([e, e], [y, y - DROP / AX_H], transform=tr, color=INK, lw=LW,
                            clip_on=False, solid_capstyle="butt")
            ax.text((XC[i] + XC[j]) / 2, y + 0.6 / AX_H, fmt_p(padj[(i, j)]), transform=tr,
                    ha="center", va="bottom", fontsize=P_PT, color=INK, clip_on=False)

        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.tick_params(axis="y", labelsize=TXT_PT, direction="out", length=TICK_LEN,
                       width=LW, pad=0.8)
        if with_legend:
            ax.set_xticks([(XC[0] + XC[1]) / 2, (XC[2] + XC[3]) / 2])
            ax.set_xticklabels(["Brain", "Liver"], fontsize=TXT_PT)
            ax.tick_params(axis="x", length=0, pad=DOT_D / 2 + 1.6)   # clear the 0 % dots
        else:
            ax.set_xticks(XC)
            ax.set_xticklabels([name for _, name, _ in CELLS], fontsize=TXT_PT, rotation=45,
                               ha="right", rotation_mode="anchor")
            ax.tick_params(axis="x", length=0, pad=DOT_D / 2 + 2.0)
        ax.set_ylabel(ylab, fontsize=TXT_PT, linespacing=1.2)
        ax.yaxis.set_label_coords(-27.4 / AX_W, 0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(LW)

        ax.text(0.5, 1 + TITLE_D / AX_H, title, transform=ax.transAxes, fontsize=TITLE_PT,
                fontweight="bold", ha="center", va="baseline")
        ax.text(-45.0 / AX_W, 1 + TITLE_D / AX_H, letter, transform=ax.transAxes,
                fontsize=LETTER_PT, fontweight="bold", ha="left", va="baseline")

        if with_legend:
            leg = ax.legend(
                handles=[Line2D([], [], marker="o", ls="", markerfacecolor=c,
                                markeredgecolor="none", markersize=DOT_D, label=lab)
                         for c, lab in ((INK, "Non-neoplastic"), (RED, "Metastatic"))],
                loc="lower left", bbox_to_anchor=(1 + 4.1 / AX_W, 0), borderaxespad=0,
                frameon=True, fancybox=False, framealpha=1, edgecolor=INK, fontsize=TXT_PT,
                handlelength=0.62, handletextpad=0.63, borderpad=0.47, labelspacing=0.77)
            leg.get_frame().set_linewidth(LW)
    return fig


os.makedirs("outputs/figures/cohort", exist_ok=True)
for stem, with_legend in (("fig_cohort_bars_pct_legend", True),
                          ("fig_cohort_bars_pct_nolegend", False)):
    fig = build(with_legend)
    fig.canvas.draw()
    k = 72.0 / fig.dpi
    for ax in fig.axes:
        bb = ax.get_window_extent()
        msg = f"    {stem}: plot {bb.width * k:.1f} x {bb.height * k:.1f} pt"
        if ax.get_legend():
            lb = ax.get_legend().get_window_extent()
            msg += f"   legend {lb.width * k:.1f} x {lb.height * k:.1f} pt (ref 61.6 x 25.3)"
        print(msg)
    for ext in ("pdf", "svg", "png"):
        q = f"outputs/figures/cohort/{stem}.{ext}"
        fig.savefig(q, dpi=600, facecolor="#ffffff", bbox_inches="tight", pad_inches=3 / 72)
        print(f"  wrote {q}")
    plt.close(fig)

with open("outputs/tables/cohort_stats_bars_pct.txt", "w") as fh:
    fh.write("Statistics printed on fig_cohort_bars_pct_{legend,nolegend} (FOV-mean endpoint).\n"
             "Kruskal-Wallis omnibus; two-sided Mann-Whitney U on all six pairs;\n"
             "Benjamini-Hochberg across those six (METHODS section 10).\n"
             "Bars = group median, whiskers = interquartile range, % of tissue area.\n\n")
    for letter, col, title, _ in PANELS:
        vals, H, p_kw, raw, padj = RESULTS[col]
        fh.write(f"=== Panel {letter}: {title.replace(SUP, '+ ')}  [{col}] ===\n")
        for (_, name, _), v in zip(CELLS, vals):
            q1, q3 = np.percentile(v, [25, 75])
            fh.write(f"  {name:21s} n={len(v):2d}  median={np.median(v):8.4f}  "
                     f"IQR {q1:.4f}-{q3:.4f}  max={v.max():.4f}\n")
        fh.write(f"  Kruskal-Wallis H={H:.3f} p={p_kw:.3g}\n")
        for (i, j), r in zip(ALL_PAIRS, raw):
            d = "drawn" if padj[(i, j)] < 0.05 else "omitted (ns)"
            fh.write(f"    {CELLS[i][1]:21s} vs {CELLS[j][1]:21s} p={r:.3g}  "
                     f"p_adj={padj[(i, j)]:.3g}  {d}\n")
        fh.write("\n")
print("  wrote outputs/tables/cohort_stats_bars_pct.txt")
