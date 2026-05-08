#!/usr/bin/env python3
"""Pairwise synteny strip plot for pangenome4 (subs + indels + rearrangements).

Same layout as synteny_pairwise.py but block widths reflect post-indel
variant segment lengths, so colinear segments may be slightly wider/narrower
than the corresponding ref segment.
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Patch

REF_SEG_LEN = 1000
N_WIN = 10
CHR_GAP = 1500
ROW_H = 0.8
PAIR_GAP = 1.2

truth = pd.read_csv("genomes4/truth.tsv", sep="\t")
variants = ["div01","div02","div05","div10","div15"]
lookup = {(r.sample, r.chr, r.slot): r for r in truth.itertuples()}

def x_for_chr(chrom):
    return 0 if chrom == "chr1" else (N_WIN*REF_SEG_LEN + CHR_GAP)

def event_color(orig_chr, orig_seg, var_chr, slot, orient):
    if orig_chr != var_chr:
        return "#3b6cb6", 0.55, "translocation"
    if orient == "-":
        return "#d23a47", 0.55, "inversion"
    if int(orig_seg) != slot:
        return "#e07a2c", 0.55, "swap"
    return "#888888", 0.30, "colinear"

fig, axes = plt.subplots(len(variants), 1, figsize=(14, 2.4*len(variants)),
                         sharex=True)

def draw_ref_strip(ax, chrom, y):
    base = x_for_chr(chrom)
    for s in range(N_WIN):
        rect = Rectangle((base + s*REF_SEG_LEN, y), REF_SEG_LEN, ROW_H,
                         facecolor="#cccccc", edgecolor="black",
                         linewidth=0.4, alpha=0.9)
        ax.add_patch(rect)

def draw_variant_strip(ax, sample, chrom, y):
    base = x_for_chr(chrom)
    for slot in range(N_WIN):
        r = lookup.get((sample, chrom, slot))
        if r is None: continue
        col, _, _ = event_color(r.orig_chr, r.orig_seg, chrom, slot, r.orient)
        hatch = "///" if r.orient == "-" else None
        rect = Rectangle((base + r.var_start, y), r.var_len, ROW_H,
                         facecolor=col, edgecolor="black",
                         linewidth=0.4, hatch=hatch, alpha=0.85)
        ax.add_patch(rect)

for vi, var in enumerate(variants):
    ax = axes[vi]
    y_top = PAIR_GAP
    y_bot = 0

    # Ref strips
    for chrom in ("chr1","chr2"):
        draw_ref_strip(ax, chrom, y_top)
        draw_variant_strip(ax, var, chrom, y_bot)

    # Sample labels
    ax.text(-700, y_top + ROW_H/2, "ref", ha="right", va="center",
            fontsize=10, fontweight="bold")
    ax.text(-700, y_bot + ROW_H/2, var, ha="right", va="center",
            fontsize=10, fontweight="bold")

    # Ribbons: each variant block to its ref origin
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
            col, alpha, _ = event_color(r.orig_chr, r.orig_seg, v_chr, slot, r.orient)
            if r.orient == "-":
                pts = [(x_r_l, y_top), (x_r_r, y_top),
                       (x_v_l, y_bot+ROW_H), (x_v_r, y_bot+ROW_H)]
            else:
                pts = [(x_r_l, y_top), (x_r_r, y_top),
                       (x_v_r, y_bot+ROW_H), (x_v_l, y_bot+ROW_H)]
            ax.add_patch(Polygon(pts, facecolor=col, alpha=alpha, edgecolor="none"))

    if vi == 0:
        ax.text(N_WIN*REF_SEG_LEN/2, y_top + ROW_H + 0.4, "chr1",
                ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.text(N_WIN*REF_SEG_LEN + CHR_GAP + N_WIN*REF_SEG_LEN/2,
                y_top + ROW_H + 0.4, "chr2",
                ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_xlim(-2000, 2*N_WIN*REF_SEG_LEN + CHR_GAP + 500)
    ax.set_ylim(-0.5, y_top + ROW_H + 1.2)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

axes[-1].set_xticks([0, N_WIN*REF_SEG_LEN, N_WIN*REF_SEG_LEN+CHR_GAP, 2*N_WIN*REF_SEG_LEN+CHR_GAP])
axes[-1].set_xticklabels(["0","10kb","0","10kb"])
axes[-1].set_xlabel("position")

legend = [
    Patch(facecolor="#cccccc", label="ref"),
    Patch(facecolor="#888888", label="colinear in variant"),
    Patch(facecolor="#e07a2c", label="within-chr swap"),
    Patch(facecolor="#d23a47", hatch="///", label="inversion"),
    Patch(facecolor="#3b6cb6", label="translocation chr1<->chr2"),
]
fig.legend(handles=legend, loc="lower center", ncol=5,
           bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=9)
fig.suptitle("Pairwise synteny vs. ref — truth (pangenome4: subs + indels + rearrangements)",
             y=0.99, fontsize=12)
plt.tight_layout()
plt.savefig("genomes4/synteny_pairwise_truth.png", dpi=140, bbox_inches="tight")
print("wrote genomes4/synteny_pairwise_truth.png")
