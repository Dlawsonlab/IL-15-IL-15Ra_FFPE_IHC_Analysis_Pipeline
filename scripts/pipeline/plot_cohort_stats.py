#!/usr/bin/env python3
"""Cohort bar charts + significance for IL-15Ra magenta area fraction, by sample type.

Writes:
- fig_magenta_bars.png : bars (group mean) + SEM + individual points + significance
  lines, on a linear % axis over the data range. Two panels: whole-slide | FOV-mean.
- fig_dab_vs_magenta.png : grouped bars, IL-15 (DAB) vs IL-15Ra (magenta) per group, on
  one shared 0-100% axis.
- outputs/tables/cohort_stats.txt : Kruskal-Wallis + pairwise Mann-Whitney U (BH-corrected).

The tests are non-parametric because the data is skewed over orders of magnitude.
"""
import math
import csv, os, math, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from itertools import combinations

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TABLE=os.path.join(ROOT,"outputs","tables","cohort_table.csv")
GROUPS=["brain","met","liver","liver_met"]   # organ-paired: brain|BCBM, then liver|liver-met
LAB={"brain":"Brain\n(control)","liver":"Normal\nliver","liver_met":"Liver\nmet","met":"BCBM"}
# One hue per ORGAN; the darker shade is the tumour/metastasis.
COL={"brain":"#9cc3ea","met":"#17529b","liver":"#f6b394","liver_met":"#b2400f"}
EDGE={"brain":"#000000","met":"#000000","liver":"#000000","liver_met":"#000000"}   # Prism: black outlines
def _nice_step(top, target_ticks=6):
    """Round tick interval (1/2/2.5/5 x 10^n) giving roughly target_ticks divisions up to top."""
    if top <= 0: return 1.0
    raw = top/float(target_ticks)
    mag = 10.0**math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m*mag: return m*mag
    return 10*mag

TICK_MAX=30      # unused; bar_panel derives the % axis top from the data
COL_DAB="#8c501e"; COL_MAG="#c8286e"   # match deconvolution display colours
FLOOR=1e-6

rows=list(csv.DictReader(open(TABLE)))
def vals(key,g):
    out=[]
    for r in rows:
        if r["control_group"]!=g: continue
        try: out.append(float(r[key]))
        except: pass
    return out
def stars(p):
    return "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"

# ---------- significance (whole-slide magenta) ----------
def sig_block(key, title):
    data=[vals(key,g) for g in GROUPS]
    lines=[f"=== {title} ({key}) ==="]
    present=[(g,d) for g,d in zip(GROUPS,data) if len(d)>=3]
    if len(present)>=2:
        H,pk=stats.kruskal(*[d for _,d in present])
        lines.append(f"Kruskal-Wallis across {len(present)} groups: H={H:.3f}, p={pk:.3g}")
    pairs=list(combinations(range(len(GROUPS)),2))
    raw=[]
    for i,j in pairs:
        a,b=data[i],data[j]
        if len(a)>=3 and len(b)>=3:
            U,p=stats.mannwhitneyu(a,b,alternative="two-sided")
            raw.append((i,j,p))
    # Benjamini-Hochberg
    m=len(raw); order=sorted(range(m),key=lambda k:raw[k][2]); adj=[None]*m
    for rank,k in enumerate(order,1):
        adj[k]=min(1.0, raw[k][2]*m/rank)
    for k in reversed(range(m-1)): adj[order[k]]=min(adj[order[k]],adj[order[k+1]])
    padj={}
    for k,(i,j,p) in enumerate(raw):
        padj[(i,j)]=adj[k]
        lines.append(f"  {GROUPS[i]:10s} vs {GROUPS[j]:10s}: p={p:.3g}  p_adj(BH)={adj[k]:.3g}  {stars(adj[k])}")
    for g,d in zip(GROUPS,data):
        if d: lines.append(f"  n({g})={len(d)}  median={st.median(d):.5g}  mean={st.mean(d):.5g}")
    return "\n".join(lines), padj

txt_ws,padj_ws = sig_block("ws_magenta_af","Whole-slide magenta area fraction")
txt_fov,padj_fov = sig_block("fov_mean_magenta_af","FOV-mean magenta area fraction")
os.makedirs(os.path.join(ROOT,"outputs","tables"),exist_ok=True)
os.makedirs(os.path.join(ROOT,"outputs","figures"),exist_ok=True)
with open(os.path.join(ROOT,"outputs","tables","cohort_stats.txt"),"w") as f:
    f.write(txt_ws+"\n\n"+txt_fov+"\n")
print(txt_ws)

# ---------- Figure 1: magenta bars + points + significance (GraphPad-Prism style) ----------
# The y axis is LINEAR in PERCENT. Bar height is then directly proportional to the value,
# which is what a reader who ignores the axis sees. The axis spans the DATA range. A 0-100%
# axis would squash all four bars into the bottom of the panel, where they read as "all
# zero". Low-signal detail belongs on a log axis, drawn as dots, not bars.
def bar_panel(ax,key,title,padj=None,logy=False,hi=None):
    ax.set_title(title,fontsize=12,fontweight="bold",color="#000")
    allpts=[]
    for gi,g in enumerate(GROUPS):
        d=[v*100 for v in vals(key,g)]            # -> percent
        if not d: continue
        if logy: d=[max(v,1e-4) for v in d]
        m=st.mean(d); sem=(st.pstdev(d)/math.sqrt(len(d))) if len(d)>1 else 0
        ax.bar(gi,m,width=0.60,color=COL[g],edgecolor="#000",linewidth=1.0,zorder=2)
        if sem>0: ax.errorbar(gi,m,yerr=sem,ecolor="#000",elinewidth=1.1,capsize=5,capthick=1.1,zorder=4)
        xs=[gi+((h*0.105)%0.28)-0.14 for h in range(len(d))]
        # white halo under each point so it stays visible on BOTH dark and light bars
        ax.scatter(xs,d,s=52,color="#ffffff",linewidths=0,zorder=5)
        ax.scatter(xs,d,s=26,color=COL[g],edgecolor="#000",linewidth=0.8,zorder=6)
        allpts+=d
    if hi is None: hi=max(allpts) if allpts else 1
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels([LAB[g]+f"\n(n={len(vals(key,g))})" for g in GROUPS],fontsize=10,color="#000")
    ax.set_ylabel("IL-15R\u03b1-positive area (% of tissue)",fontsize=10.5,color="#000")
    # Prism look: white panel, black L-shaped axes, outward ticks, NO gridlines
    ax.set_facecolor("#ffffff")
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    for sp in ("left","bottom"): ax.spines[sp].set_linewidth(1.1); ax.spines[sp].set_color("#000")
    ax.tick_params(direction="out",length=4.5,width=1.1,colors="#000",labelsize=9.5)
    ax.grid(False)
    if logy:
        ax.set_yscale("log"); ax.set_ylim(1e-4*0.8, hi*40)
    # Significance lines are flat and stack ABOVE the data, with no tall bracket legs.
    if padj:
        pairs=sorted(padj.keys(), key=lambda t:(abs(t[1]-t[0]), t[0]))
        if logy:
            y=hi*1.9
            for (i,j) in pairs:
                p=padj[(i,j)]
                ax.plot([i,j],[y,y],color="#000",lw=1.0,zorder=7,solid_capstyle="butt")
                ax.text((i+j)/2,y*1.06,stars(p),ha="center",va="bottom",fontsize=10,color="#000")
                y*=2.35
            ax.set_ylim(1e-4*0.8, y*1.4)
        else:
            # The bracket stack is compact. A small step keeps six brackets from filling the
            # panel and squashing the data into the lower half.
            step=hi*0.038; y=hi*1.03
            for (i,j) in pairs:
                p=padj[(i,j)]
                ax.plot([i,j],[y,y],color="#000",lw=1.0,zorder=7,solid_capstyle="butt")
                # The star sits ON the line with a white knockout, so the line breaks around
                # it. A label above each line would need vertical clearance per bracket and
                # would crowd the line above. This way the stack packs tightly.
                ax.text((i+j)/2,y,stars(p),ha="center",va="center",fontsize=8.5,color="#000",
                        zorder=8,bbox=dict(boxstyle="square,pad=0.18",fc="#ffffff",ec="none"))
                y+=step
            # The axis top is DERIVED from the data plus the significance stack, then rounded
            # up to a round tick. A fixed top leaves the data under a wide empty band whenever
            # exclusions lower the maximum.
            # The last labelled tick IS the axis limit, so the spine ends exactly on a labelled
            # tick. That leaves no unlabelled stub above and no truncation of the data.
            top_needed = y + step*0.30
            tick = _nice_step(top_needed)
            top = math.ceil(top_needed/tick)*tick
            ax.set_ylim(0, top)
            ax.set_yticks([i*tick for i in range(int(round(top/tick))+1)])
            ax.spines["left"].set_bounds(0, top)
    elif not logy:
        ax.set_ylim(0, hi*1.08)

def make_fig(logy,fname,note):
    fig,(a1,a2)=plt.subplots(1,2,figsize=(11,5.6),sharey=True)
    fig.patch.set_facecolor("#ffffff")
    HI=max([v*100 for k in ("ws_magenta_af","fov_mean_magenta_af") for g in GROUPS for v in vals(k,g)] or [1])
    bar_panel(a1,"ws_magenta_af","Whole-slide",padj_ws,logy=logy,hi=HI)
    bar_panel(a2,"fov_mean_magenta_af","FOV mean (per-sample mean of its FOVs)",padj_fov,logy=logy,hi=HI)
    a2.set_ylabel("")
    fig.suptitle("IL-15R\u03b1 IHC \u2014 positive area fraction by sample type",
                 fontsize=13.5,fontweight="bold",y=0.985,color="#000")
    fig.text(0.5,0.030,"Bar = group mean, error bar = SEM, points = individual samples "
             "(FOV panel: each point = that sample\u2019s mean across its own FOVs).",
             ha="center",fontsize=8.5,color="#333")
    fig.text(0.5,0.008,"Colour = organ (blue brain, orange liver); DARKER shade = tumour/metastasis. "
             "Mann-Whitney U, BH-corrected: *<.05 **<.01 ***<.001. "+note,
             ha="center",fontsize=8.5,color="#333")
    fig.tight_layout(rect=[0,0.062,1,0.955])
    fig.savefig(os.path.join(ROOT,"outputs","figures","cohort",fname),dpi=200,facecolor="#ffffff")
    print("wrote",fname)

make_fig(False,"fig_magenta_bars.png","Linear % axis over the data range.")
# Do not add a log-axis version of this bar figure (make_fig(True, ...)). A bar encodes
# magnitude by its length from a baseline, and a log axis has no meaningful baseline. The
# bars would start wherever the axis is cut, and the liver:brain difference would read as
# ~5x when it is ~4000x. Draw a log view as a dot plot instead. A dot plot has no baseline
# to distort, and it shows median/IQR to match the non-parametric tests.

# ---------- Figure 2: DAB vs magenta (GraphPad-Prism style) ----------
# This figure uses chromogen colours on purpose. Here the colour encodes the STAIN
# (brown = DAB/IL-15, magenta = Vector Red/IL-15Ra), and the x axis carries the organ.
# Do not stack the two channels. DAB-positive and magenta-positive pixels are not mutually
# exclusive, because a pixel can be both. A stacked bar would not show parts of a whole, and
# its total would be meaningless. The figure therefore uses one grouped panel.
fig2,ax=plt.subplots(figsize=(8.2,5.6))
fig2.patch.set_facecolor("#ffffff"); ax.set_facecolor("#ffffff")
w=0.34
for gi,g in enumerate(GROUPS):
    dd=[v*100 for v in vals("ws_brown_af",g)]      # DAB, percent
    dm=[v*100 for v in vals("ws_magenta_af",g)]    # magenta, percent
    for off,d,col in ((-w/2,dd,COL_DAB),(w/2,dm,COL_MAG)):
        if not d: continue
        m=st.mean(d); sem=(st.pstdev(d)/math.sqrt(len(d))) if len(d)>1 else 0
        lab = ("IL-15 (DAB)" if col==COL_DAB else "IL-15Rα (Vector Red)") if gi==0 else ""
        ax.bar(gi+off,m,w,color=col,edgecolor="#000",linewidth=1.0,zorder=2,label=lab)
        if sem>0: ax.errorbar(gi+off,m,yerr=sem,ecolor="#000",elinewidth=1.1,capsize=4,capthick=1.1,zorder=4)
        xs=[gi+off+((h*0.055)%0.15)-0.075 for h in range(len(d))]
        ax.scatter(xs,d,s=40,color="#ffffff",linewidths=0,zorder=5)          # halo
        ax.scatter(xs,d,s=18,color=col,edgecolor="#000",linewidth=0.7,zorder=6)
ax.set_ylim(0,100); ax.set_yticks(range(0,101,20))
ax.set_xticks(range(len(GROUPS)))
ax.set_xticklabels([LAB[g]+f"\n(n={len(vals('ws_brown_af',g))})" for g in GROUPS],fontsize=10,color="#000")
ax.set_ylabel("positive area (% of tissue)",fontsize=10.5,color="#000")
for sp in ("top","right"): ax.spines[sp].set_visible(False)
for sp in ("left","bottom"): ax.spines[sp].set_linewidth(1.1); ax.spines[sp].set_color("#000")
ax.tick_params(direction="out",length=4.5,width=1.1,colors="#000",labelsize=9.5)
ax.grid(False)
ax.legend(fontsize=9.5,frameon=False,loc="upper center",bbox_to_anchor=(0.5,1.045),ncol=2,handlelength=1.6,columnspacing=2.0)
fig2.suptitle("IL-15 (DAB) vs IL-15Rα (Vector Red) positive area",fontsize=13,fontweight="bold",y=0.985,color="#000")
fig2.text(0.5,0.012,
    "Bars = group mean, error bars = SEM, points = individual samples. Same 0-100% axis for both channels.\n"
    "DAB is positive across ~47-81% of tissue in EVERY group, brain control included: it saturates and does not\n"
    "discriminate. IL-15R\u03b1 is the discriminating channel.",
    ha="center", va="bottom", fontsize=8.5, color="#333", linespacing=1.6)
fig2.tight_layout(rect=[0,0.105,1,0.945])
fig2.savefig(os.path.join(ROOT,"outputs","figures","cohort","fig_dab_vs_magenta.png"),dpi=200,facecolor="#ffffff")
print("wrote fig_dab_vs_magenta.png")
