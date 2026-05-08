#!/usr/bin/env python3
"""Pairwise synteny plot for pangenome4 in Synvisio-like style:
two horizontal genome tracks per pair, smooth bezier ribbons connecting
homologous blocks. Inversions are crossed (twisted) ribbons."""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch, Patch
from matplotlib.path import Path

REF_SEG_LEN = 1000
N_WIN = 10
CHR_GAP = 1500
ROW_H = 0.6
PAIR_VSPAN = 2.4   # vertical distance between top and bottom strips

truth = pd.read_csv("genomes4/truth.tsv", sep="\t")
variants = ["div01","div02","div05","div10","div15"]
lookup = {(r.sample, r.chr, r.slot): r for r in truth.itertuples()}

CHR_COLOR = {"chr1": "#2c8cb6", "chr2": "#c66b3a"}  # ref-chromosome palettes

def x_for_chr(chrom):
    return 0 if chrom == "chr1" else (N_WIN*REF_SEG_LEN + CHR_GAP)

def event_classify(orig_chr, orig_seg, var_chr, slot, orient):
    if orig_chr != var_chr:
        return "translocation"
    if orient == "-":
        return "inversion"
    if int(orig_seg) != slot:
        return "swap"
    return "colinear"

def bezier_ribbon(x_top_l, x_top_r, y_top, x_bot_l, x_bot_r, y_bot, twist=False):
    """Build a closed bezier ribbon polygon between two horizontal segments.
    If twist=True, swap the two bottom corners (used for inversions)."""
    if twist:
        x_bot_l, x_bot_r = x_bot_r, x_bot_l
    mid = (y_top + y_bot) / 2.0
    verts = [
        (x_top_l, y_top),                                # M
        (x_top_l, mid), (x_bot_l, mid), (x_bot_l, y_bot),  # C ... to bottom-left
        (x_bot_r, y_bot),                                # L bottom-right
        (x_bot_r, mid), (x_top_r, mid), (x_top_r, y_top),  # C back to top-right
        (x_top_l, y_top),                                # close
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    return Path(verts, codes)

fig, axes = plt.subplots(len(variants), 1,
                         figsize=(14, 2.6*len(variants)), sharex=True)

def draw_track(ax, sample, chrom, y, ref_view=False):
    base = x_for_chr(chrom)
    chrom_col = CHR_COLOR[chrom]
    if ref_view:
        # ref: 10 evenly-spaced segments tinted by chromosome
        for s in range(N_WIN):
            ax.add_patch(Rectangle((base + s*REF_SEG_LEN, y), REF_SEG_LEN, ROW_H,
                                   facecolor=chrom_col, edgecolor="black",
                                   linewidth=0.4, alpha=0.85))
        return
    # Variant: blocks colored by their ORIGIN chromosome (so a chr1-origin block
    # in a chr2 row immediately reads as a translocation)
    for slot in range(N_WIN):
        r = lookup.get((sample, chrom, slot))
        if r is None: continue
        col = CHR_COLOR[r.orig_chr]
        hatch = "///" if r.orient == "-" else None
        ax.add_patch(Rectangle((base + r.var_start, y), r.var_len, ROW_H,
                               facecolor=col, edgecolor="black",
                               linewidth=0.4, hatch=hatch, alpha=0.85))

for vi, var in enumerate(variants):
    ax = axes[vi]
    y_top = PAIR_VSPAN
    y_bot = 0

    # Top = ref, bottom = variant
    for chrom in ("chr1","chr2"):
        draw_track(ax, "ref", chrom, y_top, ref_view=True)
        draw_track(ax, var, chrom, y_bot, ref_view=False)

    # Sample labels
    ax.text(-700, y_top + ROW_H/2, "ref", ha="right", va="center",
            fontsize=10, fontweight="bold")
    ax.text(-700, y_bot + ROW_H/2, var, ha="right", va="center",
            fontsize=10, fontweight="bold")

    # Bezier ribbons (variant -> ref of origin)
    for v_chr in ("chr1","chr2"):
        v_base = x_for_chr(v_chr)
        for slot in range(N_WIN):
            r = lookup.get((var, v_chr, slot))
            if r is None: continue
            r_base = x_for_chr(r.orig_chr)
            x_v_l = v_base + r.var_start
            x_v_r = x_v_l + r.var_len
            x_r_l = r_base + int(r.orig_seg)*REF_SEG_LEN
            x_r_r = x_r_l + REF_SEG_LEN
            cls = event_classify(r.orig_chr, r.orig_seg, v_chr, slot, r.orient)
            ribbon_color = CHR_COLOR[r.orig_chr]
            alpha = 0.45 if cls != "colinear" else 0.20
            twist = (r.orient == "-")
            path = bezier_ribbon(x_r_l, x_r_r, y_top,
                                 x_v_l, x_v_r, y_bot + ROW_H,
                                 twist=twist)
            ax.add_patch(PathPatch(path, facecolor=ribbon_color,
                                   edgecolor="none", alpha=alpha))

    if vi == 0:
        ax.text(N_WIN*REF_SEG_LEN/2, y_top + ROW_H + 0.3, "chr1",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=CHR_COLOR["chr1"])
        ax.text(N_WIN*REF_SEG_LEN + CHR_GAP + N_WIN*REF_SEG_LEN/2,
                y_top + ROW_H + 0.3, "chr2",
                ha="center", va="bottom", fontsize=11, fontweight="bold",
                color=CHR_COLOR["chr2"])

    ax.set_xlim(-2000, 2*N_WIN*REF_SEG_LEN + CHR_GAP + 500)
    ax.set_ylim(-0.4, y_top + ROW_H + 1.0)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

axes[-1].set_xticks([0, N_WIN*REF_SEG_LEN, N_WIN*REF_SEG_LEN+CHR_GAP, 2*N_WIN*REF_SEG_LEN+CHR_GAP])
axes[-1].set_xticklabels(["0","10kb","0","10kb"])
axes[-1].set_xlabel("position")

legend = [
    Patch(facecolor=CHR_COLOR["chr1"], label="chr1 origin"),
    Patch(facecolor=CHR_COLOR["chr2"], label="chr2 origin"),
    Patch(facecolor="grey", hatch="///", label="inversion (twisted ribbon)"),
]
fig.legend(handles=legend, loc="lower center", ncol=3,
           bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=10)
fig.suptitle("Pairwise synteny — pangenome4 (Synvisio-style ribbons)",
             y=0.99, fontsize=12)
plt.tight_layout()
plt.savefig("genomes4/synteny_pairwise_synvisio.png", dpi=140, bbox_inches="tight")
print("wrote genomes4/synteny_pairwise_synvisio.png")
