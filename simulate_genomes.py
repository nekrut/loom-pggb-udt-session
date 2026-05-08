#!/usr/bin/env python3
"""
Simulate a small pangenome for PGGB testing.

- Reference genome: 2 chromosomes, 1 Mb each (flankLen = 0).
- 5 variant genomes: each derived from the reference at substitution rates
  1%, 2%, 5%, 10%, 15%.
- Each variant mutates each chromosome independently from its corresponding
  reference chromosome, so chr1 in variants is homologous to chr1 in ref,
  and likewise for chr2.

Outputs one multi-FASTA per genome under genomes/.
"""

from random import seed as random_seed
from echydna import EchyDna
import os

CHR_LEN = 1_000_000
N_CHR = 2
DIVERGENCES = [0.01, 0.02, 0.05, 0.10, 0.15]
OUTDIR = "genomes"

random_seed(42)
EchyDna.background = "AAACCGGTTT"

os.makedirs(OUTDIR, exist_ok=True)

# Reference: 2 independent chromosomes
ref_chrs = [EchyDna(CHR_LEN) for _ in range(N_CHR)]

# Write reference as multi-FASTA
ref_path = os.path.join(OUTDIR, "genome_ref.fa")
with open(ref_path, "w") as fh:
    for i, chrom in enumerate(ref_chrs, start=1):
        chrom.fasta(fh, name=f"ref_chr{i}", wrap=100)
print(f"wrote {ref_path}  ({N_CHR} x {CHR_LEN} bp)")

# Variants
for div in DIVERGENCES:
    label = f"div{int(div*100):02d}"  # div01, div02, div05, div10, div15
    out_path = os.path.join(OUTDIR, f"genome_{label}.fa")
    with open(out_path, "w") as fh:
        for i, chrom in enumerate(ref_chrs, start=1):
            mutated = chrom.mutate(subRate=div)
            mutated.fasta(fh, name=f"{label}_chr{i}", wrap=100)
    print(f"wrote {out_path}  (subRate={div})")
