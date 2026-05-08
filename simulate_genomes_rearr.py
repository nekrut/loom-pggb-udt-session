#!/usr/bin/env python3
"""
Simulated pangenome with substitutions AND rearrangements.

- Reference: 2 chromosomes, each built from 10 segments of 1 kb (= 10 kb).
- 5 variants at substitution rates 1/2/5/10/15%; each segment mutated.
- Rearrangement events scale with divergence:
    div01: 0, div02: 1, div05: 3, div10: 6, div15: 9
- Each event is one of (1/3 each): segment swap (within-chr), inversion
  (reverse-complement one segment within-chr), or translocation
  (swap one segment between chr1 and chr2).
"""

from random import seed as random_seed, random, choice, sample, randrange
from echydna import EchyDna
import os, re

SEG_LEN  = 1_000
N_SEG    = 10                 # per chromosome
N_CHR    = 2
DIVERGENCES = [0.01, 0.02, 0.05, 0.10, 0.15]
OUTDIR   = "genomes3"

# Number of rearrangement events per variant
def n_events(div):
    return int(round(60 * div))   # div01:1? -> 0.6 -> 1 ; let's clamp:
# Actually use the discussed mapping explicitly:
EVENTS = {0.01: 0, 0.02: 1, 0.05: 3, 0.10: 6, 0.15: 9}

random_seed(42)
EchyDna.background = "AAACCGGTTT"
os.makedirs(OUTDIR, exist_ok=True)

# Reference: list of N_CHR chromosomes, each a list of N_SEG segments
ref = [[EchyDna(SEG_LEN) for _ in range(N_SEG)] for _ in range(N_CHR)]

def write_genome(path, chrs, sample_name):
    with open(path, "w") as fh:
        for i, segs in enumerate(chrs, start=1):
            # concatenate all segments into one sequence
            seq = segs[0]
            for s in segs[1:]:
                seq = seq + s
            seq.fasta(fh, name=f"{sample_name}_chr{i}", wrap=100)

# Reference output
write_genome(f"{OUTDIR}/genome_ref.fa", ref, "ref")
print(f"wrote ref ({N_CHR} chr x {N_SEG*SEG_LEN} bp)")

# Variants
for div in DIVERGENCES:
    label = f"div{int(div*100):02d}"
    # Mutate each segment at this rate (independent per segment)
    var = [[seg.mutate(subRate=div) for seg in chrom] for chrom in ref]

    n = EVENTS[div]
    counts = {"swap":0, "inv":0, "tloc":0}
    for _ in range(n):
        op = choice(["swap","inv","tloc"])
        if op == "swap":
            c = randrange(N_CHR)
            i, j = sample(range(N_SEG), 2)
            var[c][i], var[c][j] = var[c][j], var[c][i]
        elif op == "inv":
            c = randrange(N_CHR)
            i = randrange(N_SEG)
            var[c][i] = -var[c][i]   # reverse-complement
        else:  # tloc: swap one segment between chr1 and chr2
            i = randrange(N_SEG)
            j = randrange(N_SEG)
            var[0][i], var[1][j] = var[1][j], var[0][i]
        counts[op] += 1

    write_genome(f"{OUTDIR}/genome_{label}.fa", var, label)
    print(f"wrote {label}  subRate={div}  events: {counts}")

# Combine into PanSN multifasta for PGGB
samples = ["ref","div01","div02","div05","div10","div15"]
with open(f"{OUTDIR}/pangenome3.fa","w") as out:
    for s in samples:
        with open(f"{OUTDIR}/genome_{s}.fa") as fh:
            for L in fh:
                if L.startswith(">"):
                    m = re.match(r">(\S+)_(chr\d+)", L)
                    out.write(f">{m.group(1)}#1#{m.group(2)}\n")
                else:
                    out.write(L)
print("wrote pangenome3.fa")
