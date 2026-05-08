#!/usr/bin/env python3
"""
Re-run the same simulation as simulate_genomes_rearr.py (same seed)
but ALSO record, for every (sample, chromosome, window_index), which
reference (chr, segment_idx) and orientation it derives from. Saves a
truth.tsv we can use for ground-truth identity.

We also record the actual mutated segment sequence so we can compute
exact Hamming identity to the un-mutated ref segment.
"""
from random import seed as random_seed, choice, sample, randrange
from echydna import EchyDna
import os, json

SEG_LEN = 1_000
N_SEG   = 10
N_CHR   = 2
DIVERGENCES = [0.01, 0.02, 0.05, 0.10, 0.15]
EVENTS = {0.01: 0, 0.02: 1, 0.05: 3, 0.10: 6, 0.15: 9}
OUTDIR = "genomes3"

random_seed(42)
EchyDna.background = "AAACCGGTTT"

# Reference: 2 chromosomes, each 10 segments
ref = [[EchyDna(SEG_LEN) for _ in range(N_SEG)] for _ in range(N_CHR)]
ref_seqs = [[str(seg.seq) for seg in chrom] for chrom in ref]

records = []  # (sample, chrom_idx, window_idx, ref_chr, ref_seg, orient, hamming_identity)

# Reference itself
for c in range(N_CHR):
    for w in range(N_SEG):
        records.append({
            "sample":"ref","chr":f"chr{c+1}","window":w,
            "ref_chr":f"chr{c+1}","ref_seg":w,"orient":"+",
            "identity":1.0
        })

def hamming_id(a,b):
    return sum(x==y for x,y in zip(a,b)) / max(len(a),1)

def revcomp(s):
    comp={"A":"T","T":"A","G":"C","C":"G","N":"N"}
    return "".join(comp.get(b,"N") for b in reversed(s))

for div in DIVERGENCES:
    label=f"div{int(div*100):02d}"
    # Each variant segment knows its (orig_chr, orig_seg, orient='+')
    var = [
        [
            {"seg": ref[c][s].mutate(subRate=div),
             "orig_chr": c, "orig_seg": s, "orient": "+"}
            for s in range(N_SEG)
        ]
        for c in range(N_CHR)
    ]
    # Apply rearrangements (track structure changes only; sequences move with them)
    n = EVENTS[div]
    for _ in range(n):
        op = choice(["swap","inv","tloc"])
        if op=="swap":
            c=randrange(N_CHR); i,j=sample(range(N_SEG),2)
            var[c][i],var[c][j] = var[c][j],var[c][i]
        elif op=="inv":
            c=randrange(N_CHR); i=randrange(N_SEG)
            seg = var[c][i]
            seg["seg"] = -seg["seg"]
            seg["orient"] = "-" if seg["orient"]=="+" else "+"
        else:
            i=randrange(N_SEG); j=randrange(N_SEG)
            var[0][i],var[1][j] = var[1][j],var[0][i]

    # Now compute truth identity for each window
    for c in range(N_CHR):
        for w in range(N_SEG):
            seg_info = var[c][w]
            ref_seq = ref_seqs[seg_info["orig_chr"]][seg_info["orig_seg"]]
            mutated_seq = str(seg_info["seg"].seq)
            # If inverted, the segment's sequence is the RC of original-mutation;
            # to compare against unmutated ref, we should compare in the
            # orientation in which mutation occurred (the un-RC'd form).
            if seg_info["orient"] == "-":
                mutated_seq = revcomp(mutated_seq)
            ident = hamming_id(ref_seq, mutated_seq)
            records.append({
                "sample":label,"chr":f"chr{c+1}","window":w,
                "ref_chr":f"chr{seg_info['orig_chr']+1}",
                "ref_seg":seg_info["orig_seg"],
                "orient":seg_info["orient"],
                "identity":ident
            })

# Write TSV
with open(f"{OUTDIR}/truth.tsv","w") as fh:
    fh.write("sample\tchr\twindow\tref_chr\tref_seg\torient\tidentity\n")
    for r in records:
        fh.write(f"{r['sample']}\t{r['chr']}\t{r['window']}\t{r['ref_chr']}\t{r['ref_seg']}\t{r['orient']}\t{r['identity']:.4f}\n")
print(f"wrote {len(records)} rows to {OUTDIR}/truth.tsv")
