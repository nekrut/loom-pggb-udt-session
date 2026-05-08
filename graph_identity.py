#!/usr/bin/env python3
"""
For each variant path, slide a 1-kb window along its sequence and compute
the fraction of bases that lie on graph nodes also visited by the
*corresponding* reference chromosome path (chr1 vs ref#1#chr1, etc.).

This is the 'graph-reconstructed identity' to compare against the truth.
"""
import os, csv

GFA = "genomes3/pangenome3.gfa"
WINDOW = 1000
N_WINDOWS = 10
OUT = "genomes3/graph_identity.tsv"

# Parse S-lines: node_id -> length
node_len = {}
paths = {}  # name -> list of (node_id, orient) tuples

with open(GFA) as fh:
    for line in fh:
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            node_len[f[1]] = len(f[2])
        elif line.startswith("P\t"):
            f = line.rstrip("\n").split("\t")
            name = f[1]
            # f[2] is comma-separated like "1+,87+,5-,..."
            seg_str = f[2].split(",")
            nodes = [(s[:-1], s[-1]) for s in seg_str]
            paths[name] = nodes

print(f"loaded {len(paths)} paths, {len(node_len)} segments")

# Reference node sets, keyed by chr suffix
ref_node_sets = {
    "chr1": set(n for n,_ in paths["ref#1#chr1"]),
    "chr2": set(n for n,_ in paths["ref#1#chr2"]),
}

samples = ["ref","div01","div02","div05","div10","div15"]

rows = []
for s in samples:
    for c in ("chr1","chr2"):
        path_name = f"{s}#1#{c}"
        if path_name not in paths:
            continue
        ref_set = ref_node_sets[c]
        # For each window of WINDOW bases along the path, accumulate
        # base counts: shared / total
        # Walk the node list keeping running offset; partition into windows.
        win_total = [0]*N_WINDOWS
        win_shared = [0]*N_WINDOWS
        offset = 0
        for nid, _ori in paths[path_name]:
            L = node_len[nid]
            shared = nid in ref_set
            # node spans [offset, offset+L)
            start = offset
            end   = offset + L
            # distribute across windows
            w0 = start // WINDOW
            w1 = (end - 1) // WINDOW
            for w in range(w0, w1+1):
                if w >= N_WINDOWS: break
                wstart = w*WINDOW
                wend   = (w+1)*WINDOW
                lo = max(start, wstart)
                hi = min(end, wend)
                bases = hi - lo
                win_total[w] += bases
                if shared:
                    win_shared[w] += bases
            offset = end
        for w in range(N_WINDOWS):
            tot = win_total[w]
            if tot == 0:
                ident = float("nan")
            else:
                ident = win_shared[w] / tot
            rows.append((s, c, w, ident, tot))

with open(OUT, "w") as fh:
    fh.write("sample\tchr\twindow\tidentity\twindow_bases\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]}\n")
print(f"wrote {len(rows)} rows to {OUT}")
