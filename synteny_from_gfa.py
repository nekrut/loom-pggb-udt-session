#!/usr/bin/env python3
"""
Reconstruct synvisio-style synteny plot from a PGGB GFA.

For each variant path (PanSN-named sample#hap#chr), split into N_CHUNKS
equal-base windows along the path. For each window, find which ref chunk
shares the most graph nodes -> infer (origin_chr, origin_chunk, orient).
Then draw the same kind of ribbon plot we made from truth.

Usage:
  python synteny_from_gfa.py GFA [REF_NAME] [N_CHUNKS] [OUT_PNG]
Defaults: REF_NAME=ref, N_CHUNKS=10, OUT_PNG=synteny_from_gfa.png
"""
import sys, re, math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch, Patch
from matplotlib.path import Path

GFA      = sys.argv[1]
REF_NAME = sys.argv[2] if len(sys.argv) > 2 else "ref"
N_CHUNKS = int(sys.argv[3]) if len(sys.argv) > 3 else 10
OUT_PNG  = sys.argv[4] if len(sys.argv) > 4 else "synteny_from_gfa.png"

# ---- 1. Parse GFA ----
node_len = {}
paths = {}   # name -> list of (node_id, orient)

with open(GFA) as fh:
    for line in fh:
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            node_len[f[1]] = len(f[2])
        elif line.startswith("P\t"):
            f = line.rstrip("\n").split("\t")
            seg = f[2].split(",")
            paths[f[1]] = [(s[:-1], s[-1]) for s in seg]

print(f"loaded {len(paths)} paths, {len(node_len)} nodes from {GFA}", file=sys.stderr)

# Group paths by sample. PanSN expected: sample#hap#chr
def parse_panSN(name):
    parts = name.split("#")
    if len(parts) >= 3:
        return parts[0], parts[-1]   # (sample, chr)
    return None, None

samples_chr = defaultdict(dict)  # sample -> chr -> path_name
for pname in paths:
    s, c = parse_panSN(pname)
    if s and c:
        samples_chr[s][c] = pname

if REF_NAME not in samples_chr:
    raise SystemExit(f"reference sample '{REF_NAME}' not found among paths: "
                     f"{sorted(samples_chr.keys())}")

chroms = sorted(samples_chr[REF_NAME].keys())
non_ref = [s for s in samples_chr if s != REF_NAME]

# ---- 2. Compute path lengths and ref chunk boundaries ----
def path_len(pname):
    return sum(node_len[n] for n,_ in paths[pname])

ref_path_len = {c: path_len(samples_chr[REF_NAME][c]) for c in chroms}

# For each ref chr, build list of (chunk_idx, set_of_nodes, set_of_(node,orient))
def build_chunks(pname, n_chunks):
    """Walk path and split into n_chunks equal-base spans. Return list of
    dicts: {idx, start, end, nodes:set, oriented:set}."""
    L = path_len(pname)
    boundaries = [int(round(i * L / n_chunks)) for i in range(n_chunks+1)]
    chunks = [{"idx": i, "start": boundaries[i], "end": boundaries[i+1],
               "nodes": set(), "oriented": set()} for i in range(n_chunks)]
    offset = 0
    for nid, ori in paths[pname]:
        L_n = node_len[nid]
        a, b = offset, offset + L_n
        # node spans [a, b); assign to all chunks it overlaps
        for ch in chunks:
            lo = max(a, ch["start"])
            hi = min(b, ch["end"])
            if hi > lo:
                ch["nodes"].add(nid)
                ch["oriented"].add((nid, ori))
        offset = b
    return chunks

ref_chunks = {c: build_chunks(samples_chr[REF_NAME][c], N_CHUNKS) for c in chroms}

# ---- 3. For each variant path, classify each chunk ----
def classify(var_chunk, var_chr):
    """Return (origin_chr, origin_chunk_idx, orient_flag)."""
    best = (None, None, 0)  # (chr, chunk_idx, overlap_count)
    for cchr in chroms:
        for ch in ref_chunks[cchr]:
            overlap = len(var_chunk["nodes"] & ch["nodes"])
            if overlap > best[2]:
                best = (cchr, ch["idx"], overlap)
    if best[0] is None:
        return None, None, "+"
    # Determine orientation: compare oriented node sets
    ref_or = ref_chunks[best[0]][best[1]]["oriented"]
    var_or = var_chunk["oriented"]
    same = sum(1 for x in var_or if x in ref_or)
    flipped = sum(1 for (n,o) in var_or
                  if (n, "+" if o == "-" else "-") in ref_or)
    orient = "-" if flipped > same else "+"
    return best[0], best[1], orient

# Build classification table: dict (sample, chr, chunk_idx) -> (orig_chr, orig_seg, orient, var_start, var_len)
results = {}
for s in [REF_NAME] + non_ref:
    for c in chroms:
        if c not in samples_chr[s]: continue
        chunks = build_chunks(samples_chr[s][c], N_CHUNKS)
        for ch in chunks:
            if s == REF_NAME:
                results[(s,c,ch["idx"])] = (c, ch["idx"], "+",
                                            ch["start"], ch["end"]-ch["start"])
            else:
                oc, os_, orient = classify(ch, c)
                results[(s,c,ch["idx"])] = (oc, os_, orient,
                                            ch["start"], ch["end"]-ch["start"])

# ---- 4. Plot (synvisio-style) ----
# Color BY EVENT TYPE (matches the truth plot)
EVENT_COLOR = {
    "colinear":      "#888888",   # grey
    "swap":          "#e07a2c",   # orange
    "inversion":     "#d23a47",   # red
    "translocation": "#3b6cb6",   # blue
    "unknown":       "#cccccc",   # light grey (no overlap)
}
CHR_COLOR = {chroms[0]: "#2c8cb6", chroms[1]: "#c66b3a"} if len(chroms) >= 2 \
            else {chroms[0]: "#2c8cb6"}

def classify_event(orig_chr, orig_seg, var_chr, slot, orient):
    if orig_chr is None: return "unknown"
    if orig_chr != var_chr: return "translocation"
    if orient == "-":      return "inversion"
    if int(orig_seg) != slot: return "swap"
    return "colinear"
CHR_GAP = max(ref_path_len.values()) // 6 + 500
ROW_H = 0.6
PAIR_VSPAN = 2.4

def x_for_chr(chrom):
    if chrom == chroms[0]: return 0
    return ref_path_len[chroms[0]] + CHR_GAP

def bezier_ribbon(x_top_l, x_top_r, y_top, x_bot_l, x_bot_r, y_bot, twist=False):
    if twist:
        x_bot_l, x_bot_r = x_bot_r, x_bot_l
    mid = (y_top + y_bot) / 2
    verts = [
        (x_top_l, y_top),
        (x_top_l, mid), (x_bot_l, mid), (x_bot_l, y_bot),
        (x_bot_r, y_bot),
        (x_bot_r, mid), (x_top_r, mid), (x_top_r, y_top),
        (x_top_l, y_top),
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    return Path(verts, codes)

fig, axes = plt.subplots(len(non_ref), 1,
                         figsize=(14, 2.6*max(1,len(non_ref))), sharex=True)
if len(non_ref) == 1: axes = [axes]

def draw_track(ax, sample, chrom, y, ref_view=False):
    base = x_for_chr(chrom)
    if ref_view:
        for i in range(N_CHUNKS):
            r = results[(sample, chrom, i)]
            ax.add_patch(Rectangle((base + r[3], y), r[4], ROW_H,
                                   facecolor=CHR_COLOR[chrom], edgecolor="black",
                                   linewidth=0.4, alpha=0.85))
        return
    for i in range(N_CHUNKS):
        r = results.get((sample, chrom, i))
        if r is None: continue
        oc, os_, orient, vstart, vlen = r
        ev = classify_event(oc, os_, chrom, i, orient)
        col = EVENT_COLOR[ev]
        hatch = "///" if ev == "inversion" else None
        ax.add_patch(Rectangle((base + vstart, y), vlen, ROW_H,
                               facecolor=col, edgecolor="black",
                               linewidth=0.4, hatch=hatch, alpha=0.85))

for vi, var in enumerate(non_ref):
    ax = axes[vi]
    y_top, y_bot = PAIR_VSPAN, 0
    for chrom in chroms:
        draw_track(ax, REF_NAME, chrom, y_top, ref_view=True)
        draw_track(ax, var, chrom, y_bot, ref_view=False)
    ax.text(-700, y_top + ROW_H/2, REF_NAME, ha="right", va="center",
            fontsize=10, fontweight="bold")
    ax.text(-700, y_bot + ROW_H/2, var, ha="right", va="center",
            fontsize=10, fontweight="bold")

    # Ribbons
    for v_chr in chroms:
        v_base = x_for_chr(v_chr)
        for i in range(N_CHUNKS):
            r = results.get((var, v_chr, i))
            if r is None: continue
            oc, os_, orient, vstart, vlen = r
            if oc is None: continue
            r_base = x_for_chr(oc)
            ref_r = results[(REF_NAME, oc, os_)]
            x_v_l = v_base + vstart;          x_v_r = x_v_l + vlen
            x_r_l = r_base + ref_r[3];         x_r_r = x_r_l + ref_r[4]
            ev = classify_event(oc, os_, v_chr, i, orient)
            twist = (ev == "inversion")
            alpha = 0.20 if ev == "colinear" else 0.55
            path = bezier_ribbon(x_r_l, x_r_r, y_top,
                                 x_v_l, x_v_r, y_bot + ROW_H, twist=twist)
            ax.add_patch(PathPatch(path, facecolor=EVENT_COLOR[ev],
                                   edgecolor="none", alpha=alpha))

    if vi == 0:
        for chrom in chroms:
            ax.text(x_for_chr(chrom) + ref_path_len[chrom]/2,
                    y_top + ROW_H + 0.3, chrom,
                    ha="center", va="bottom", fontsize=11, fontweight="bold",
                    color=CHR_COLOR[chrom])

    max_x = max(x_for_chr(c) + ref_path_len[c] for c in chroms)
    ax.set_xlim(-2000, max_x + 500)
    ax.set_ylim(-0.4, y_top + ROW_H + 1.0)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

axes[-1].set_xticks([0, ref_path_len[chroms[0]],
                     x_for_chr(chroms[1]) if len(chroms)>1 else 0,
                     max(x_for_chr(c)+ref_path_len[c] for c in chroms)])
axes[-1].set_xticklabels([f"0",f"{ref_path_len[chroms[0]]/1000:.0f}kb",
                          "0", f"{ref_path_len[chroms[-1]]/1000:.0f}kb"]
                         if len(chroms)>1 else
                         ["0", f"{ref_path_len[chroms[0]]/1000:.0f}kb"])
axes[-1].set_xlabel("position")

legend = [
    Patch(facecolor=EVENT_COLOR["colinear"],      label="colinear"),
    Patch(facecolor=EVENT_COLOR["swap"],          label="within-chr swap"),
    Patch(facecolor=EVENT_COLOR["inversion"], hatch="///", label="inversion"),
    Patch(facecolor=EVENT_COLOR["translocation"], label="translocation chr1<->chr2"),
    Patch(facecolor=EVENT_COLOR["unknown"],       label="unrecovered (no shared nodes)"),
]
fig.legend(handles=legend, loc="lower center", ncol=len(legend),
           bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=10)
fig.suptitle(f"Pairwise synteny reconstructed from PGGB graph "
             f"(N={N_CHUNKS} chunks/chr; ref='{REF_NAME}')",
             y=0.99, fontsize=12)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=140, bbox_inches="tight")
print(f"wrote {OUT_PNG}", file=sys.stderr)
