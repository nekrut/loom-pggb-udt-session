#!/usr/bin/env python3
"""
Pairwise synteny strip plot: each variant vs. ref.

Layout (riparian-style):
  For each variant in [div01, div02, div05, div10, div15]:
    A pair of horizontal strips:
      - top    : REF chr1 (left) + ref chr2 (right) with small gap
      - bottom : VARIANT chr1 (left) + variant chr2 (right) with same gap
    Ribbons connect each variant block to its origin in ref:
      - colinear (same chr, same position, same orient) = grey
      - inversion                                       = red
      - translocation (chr1 <-> chr2)                   = blue
      - within-chr swap                                 = orange
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

WIN, N_WIN = 1000, 10
CHR_GAP = 1500   # gap between chr1 and chr2 panels on a row
ROW_H = 0.8
PAIR_GAP = 1.2   # gap between top (ref) and bottom (variant) within a pair
PAIR_BLOCK = 4.5 # vertical space per variant comparison

truth = pd.read_csv("genomes3/truth.tsv", sep="\t")
variants = ["div01","div02","div05","div10","div15"]

lookup = {(r.sample, r.chr, r.window): r for r in truth.itertuples()}

def x_for(chrom, w):
    """X coordinate of the LEFT edge of window w in chromosome chrom."""
    base = 0 if chrom == "chr1" else (N_WIN*WIN + CHR_GAP)
    return base + w*WIN

def event_color(r_top_chr, r_top_seg, r_top_orient,
                v_chr, v_w, v_orient):
    """Color the ribbon by event type."""
    if r_top_chr != v_chr:
        return ("#3b6cb6", 0.55, "translocation")   # blue
    if v_orient == "-":
        return ("#d23a47", 0.55, "inversion")       # red
    if int(r_top_seg) != v_w:
        return ("#e07a2c", 0.55, "swap")            # orange
    return ("#888888", 0.30, "colinear")            # grey

fig, axes = plt.subplots(len(variants), 1, figsize=(14, 2.4*len(variants)),
                         sharex=True)

# Helper to draw a strip of N_WIN blocks at vertical y, for a given sample/chr
def draw_strip(ax, sample, chrom, y, label_left=None, label_color="black"):
    for w in range(N_WIN):
        r = lookup.get((sample, chrom, w))
        if r is None: continue
        x = x_for(chrom, w)
        # block fill: light grey if colinear, otherwise color by event
        if sample == "ref":
            fc = "#cccccc"
            hatch = None
        else:
            col, _, _ = event_color(r.ref_chr, r.ref_seg, "+",
                                    chrom, w, r.orient)
            fc = col
            hatch = "///" if r.orient == "-" else None
        rect = Rectangle((x, y), WIN, ROW_H,
                         facecolor=fc, edgecolor="black",
                         linewidth=0.4, hatch=hatch, alpha=0.85)
        ax.add_patch(rect)
    if label_left is not None:
        # left-most position depends on chr1 panel start (x=0)
        ax.text(-700, y + ROW_H/2, label_left, ha="right", va="center",
                fontsize=10, fontweight="bold", color=label_color)

for vi, var in enumerate(variants):
    ax = axes[vi]
    y_top = PAIR_GAP        # ref strip
    y_bot = 0               # variant strip

    # Draw chromosome dividers / labels
    for chrom in ("chr1","chr2"):
        # ref strip
        draw_strip(ax, "ref", chrom, y_top,
                   label_left=("ref" if chrom=="chr1" else None))
        # variant strip
        draw_strip(ax, var, chrom, y_bot,
                   label_left=(var if chrom=="chr1" else None))

    # Draw ribbons: for each variant block, draw polygon to its origin in ref
    for v_chr in ("chr1","chr2"):
        for v_w in range(N_WIN):
            r = lookup.get((var, v_chr, v_w))
            if r is None: continue
            x_v_l = x_for(v_chr, v_w);     x_v_r = x_v_l + WIN
            x_r_l = x_for(r.ref_chr, int(r.ref_seg)); x_r_r = x_r_l + WIN
            col, alpha, _ = event_color(r.ref_chr, r.ref_seg, "+",
                                        v_chr, v_w, r.orient)
            # if inversion, cross the ribbon (swap top corners)
            if r.orient == "-":
                pts = [(x_r_l, y_top), (x_r_r, y_top),
                       (x_v_l, y_bot+ROW_H), (x_v_r, y_bot+ROW_H)]
            else:
                pts = [(x_r_l, y_top), (x_r_r, y_top),
                       (x_v_r, y_bot+ROW_H), (x_v_l, y_bot+ROW_H)]
            poly = Polygon(pts, facecolor=col, alpha=alpha,
                           edgecolor="none")
            ax.add_patch(poly)

    # Chromosome labels (above ref strip on first row only)
    if vi == 0:
        ax.text(N_WIN*WIN/2, y_top + ROW_H + 0.4, "chr1",
                ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(N_WIN*WIN + CHR_GAP + N_WIN*WIN/2, y_top + ROW_H + 0.4, "chr2",
                ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_xlim(-2000, 2*N_WIN*WIN + CHR_GAP + 500)
    ax.set_ylim(-0.5, y_top + ROW_H + 1.2)
    ax.set_yticks([]); ax.set_xticks([])
    for s in ax.spines.values(): s.set_visible(False)

axes[-1].set_xticks([0, N_WIN*WIN, N_WIN*WIN+CHR_GAP, 2*N_WIN*WIN+CHR_GAP])
axes[-1].set_xticklabels(["0","10kb","0","10kb"])
axes[-1].set_xlabel("position")

# Legend
from matplotlib.patches import Patch
legend = [
    Patch(facecolor="#cccccc", label="ref / colinear (grey)"),
    Patch(facecolor="#888888", label="colinear in variant"),
    Patch(facecolor="#e07a2c", label="within-chr swap"),
    Patch(facecolor="#d23a47", hatch="///", label="inversion"),
    Patch(facecolor="#3b6cb6", label="translocation chr1 <-> chr2"),
]
fig.legend(handles=legend, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.02),
           frameon=False, fontsize=9)
fig.suptitle("Pairwise synteny vs. ref — truth (each variant compared to ref)",
             y=0.99, fontsize=12)
plt.tight_layout()
plt.savefig("genomes3/synteny_pairwise_truth.png", dpi=140, bbox_inches="tight")
print("wrote genomes3/synteny_pairwise_truth.png")
