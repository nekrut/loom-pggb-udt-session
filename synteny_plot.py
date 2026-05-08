#!/usr/bin/env python3
"""
Synteny strip plot for the simulated pangenome (truth version).

- Two panels: chr1 (top), chr2 (bottom)
- Each sample is one horizontal strip of 10 colored 1-kb blocks
- Color encodes the (ref_chr, ref_seg) origin of each block
  (so chr1-origin blocks living on a chr2 strip = translocations)
- Inversions are hatched
- Ribbons between adjacent sample rows connect blocks of matching origin
  (straight = same x range, crossing = swap, ribbons jumping between panels = translocation)
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors

WIN = 1000
N_WIN = 10
truth = pd.read_csv("genomes3/truth.tsv", sep="\t")
samples = ["ref","div01","div02","div05","div10","div15"]

# Color: 20 colors total (chr1 segs = warm/red palette, chr2 segs = cool/blue palette)
# Use distinct hues so it's easy to follow ribbons
warm = plt.cm.YlOrRd(np.linspace(0.3, 0.95, N_WIN)) if False else plt.cm.tab20(range(20))
# Build palette: chr1 seg 0-9 -> first 10 of tab20, chr2 seg 0-9 -> last 10
import numpy as np
palette = {}
tab20 = plt.cm.tab20.colors
for s in range(N_WIN):
    palette[("chr1", s)] = tab20[s]            # 0..9
    palette[("chr2", s)] = tab20[10 + s]       # 10..19

def block_color(ref_chr, ref_seg):
    return palette[(ref_chr, int(ref_seg))]

# Build geometry
# panel y-coords: chr1 in top axes, chr2 in bottom axes
# samples within a panel: y from top to bottom
ROW_H = 0.7
ROW_GAP = 0.5  # vertical gap between sample rows

fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
fig.subplots_adjust(hspace=0.15)

def y_for_sample(sample_idx):
    # row 0 at top
    return -(sample_idx * (ROW_H + ROW_GAP))

# Lookup: (sample, chr, window) -> row dict
lookup = {(r.sample, r.chr, r.window): r for r in truth.itertuples()}

for panel_idx, chrom in enumerate(["chr1","chr2"]):
    ax = axes[panel_idx]
    # Draw blocks for each sample in this chromosome panel
    for si, s in enumerate(samples):
        y = y_for_sample(si)
        # sample label
        ax.text(-200, y + ROW_H/2, s, ha="right", va="center",
                fontsize=10, fontweight="bold")
        for w in range(N_WIN):
            r = lookup.get((s, chrom, w))
            if r is None: continue
            x = w * WIN
            color = block_color(r.ref_chr, r.ref_seg)
            hatch = "///" if r.orient == "-" else None
            rect = Rectangle((x, y), WIN, ROW_H,
                             facecolor=color, edgecolor="black",
                             linewidth=0.5, hatch=hatch)
            ax.add_patch(rect)
            # tiny label inside block: origin
            label = f"{r.ref_chr[-1]}{r.ref_seg}{'-' if r.orient=='-' else ''}"
            ax.text(x + WIN/2, y + ROW_H/2, label,
                    ha="center", va="center", fontsize=7, color="black")
    # Draw ribbons between adjacent rows for blocks of matching origin
    # Ribbon connects bottom edge of block in row i to top edge of block in row i+1
    for si in range(len(samples)-1):
        s_top = samples[si]
        s_bot = samples[si+1]
        y_top = y_for_sample(si)            # bottom of top row = y_top
        y_bot = y_for_sample(si+1) + ROW_H  # top of bottom row
        # For each block in s_top, find matching block in s_bot (same chromosome panel only)
        for w_top in range(N_WIN):
            r_top = lookup.get((s_top, chrom, w_top))
            if r_top is None: continue
            # find w in s_bot where origin matches
            for w_bot in range(N_WIN):
                r_bot = lookup.get((s_bot, chrom, w_bot))
                if r_bot is None: continue
                if (r_bot.ref_chr == r_top.ref_chr and
                    int(r_bot.ref_seg) == int(r_top.ref_seg)):
                    # draw a trapezoid
                    x_top_l = w_top*WIN; x_top_r = (w_top+1)*WIN
                    x_bot_l = w_bot*WIN; x_bot_r = (w_bot+1)*WIN
                    poly = Polygon(
                        [(x_top_l, y_top), (x_top_r, y_top),
                         (x_bot_r, y_bot), (x_bot_l, y_bot)],
                        facecolor=block_color(r_top.ref_chr, r_top.ref_seg),
                        alpha=0.20, edgecolor="none")
                    ax.add_patch(poly)

    ax.set_xlim(-1500, N_WIN*WIN + 200)
    ax.set_ylim(y_for_sample(len(samples)-1) - 0.2, ROW_H + 0.2)
    ax.set_yticks([])
    ax.set_title(f"{chrom}", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("position (bp)" if panel_idx == 1 else "")
    for spine in ("top","right","left"):
        ax.spines[spine].set_visible(False)

# Cross-panel translocation ribbons: connect adjacent rows where a block's
# origin chromosome differs from the panel's chromosome — but we already show
# this via color + label inside each block. To keep the figure readable we
# don't draw cross-panel polygons (would be a tangle).

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=palette[("chr1", s)], label=f"chr1 seg{s}") for s in range(N_WIN)] + \
                  [Patch(facecolor=palette[("chr2", s)], label=f"chr2 seg{s}") for s in range(N_WIN)]
fig.legend(handles=legend_elements, loc="center right",
           bbox_to_anchor=(1.05, 0.5), ncol=1, fontsize=8, frameon=False,
           title="Block origin\n(in ref)")

fig.suptitle("Synteny — truth (block color = origin in ref; hatched = inverted)",
             fontsize=12, y=0.95)
plt.savefig("genomes3/synteny_truth.png", dpi=140, bbox_inches="tight")
print("wrote genomes3/synteny_truth.png")
