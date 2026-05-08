#!/usr/bin/env python3
"""
Same as simulate_genomes_rearr.py but using echydna's built-in indel support.

- subRate = divergence (1/2/5/10/15%)
- indelRate = subRate / 10
- indels = 5  (uniform -5..+5, excluding 0)
- Rearrangements: same schedule (swap / inversion / translocation)
"""
from random import seed as random_seed, choice, sample, randrange
from echydna import EchyDna
import os, re

SEG_LEN = 100_000
N_SEG   = 10
N_CHR   = 2
DIVERGENCES = [0.01, 0.02, 0.05, 0.10, 0.15]
EVENTS = {0.01: 0, 0.02: 1, 0.05: 3, 0.10: 6, 0.15: 9}
INDEL_RATIO = 0.10
INDEL_MAX_LEN = 5
OUTDIR = "genomes5"

random_seed(42)
EchyDna.background = "AAACCGGTTT"
os.makedirs(OUTDIR, exist_ok=True)

ref = [[EchyDna(SEG_LEN) for _ in range(N_SEG)] for _ in range(N_CHR)]

def write_genome(path, chrs, sample_name):
    with open(path, "w") as fh:
        for i, segs in enumerate(chrs, start=1):
            seq = segs[0]
            for s in segs[1:]:
                seq = seq + s
            seq.fasta(fh, name=f"{sample_name}_chr{i}", wrap=100)

write_genome(f"{OUTDIR}/genome_ref.fa", ref, "ref")
print("wrote ref")

provenance = []  # for record-keeping

for div in DIVERGENCES:
    label = f"div{int(div*100):02d}"
    indel_rate = div * INDEL_RATIO
    # Mutate each ref segment with subs + indels
    var = [
        [
            ref[c][s].mutate(subRate=div, indelRate=indel_rate, indels=INDEL_MAX_LEN)
            for s in range(N_SEG)
        ]
        for c in range(N_CHR)
    ]
    # Track origin metadata so we can keep producing a synteny-truth file
    meta = [
        [{"orig_chr": c, "orig_seg": s, "orient": "+",
          "len": len(var[c][s])} for s in range(N_SEG)]
        for c in range(N_CHR)
    ]

    n = EVENTS[div]
    counts = {"swap":0,"inv":0,"tloc":0}
    for _ in range(n):
        op = choice(["swap","inv","tloc"])
        if op == "swap":
            c = randrange(N_CHR); i,j = sample(range(N_SEG), 2)
            var[c][i], var[c][j] = var[c][j], var[c][i]
            meta[c][i], meta[c][j] = meta[c][j], meta[c][i]
        elif op == "inv":
            c = randrange(N_CHR); i = randrange(N_SEG)
            var[c][i] = -var[c][i]
            meta[c][i]["orient"] = "-" if meta[c][i]["orient"] == "+" else "+"
        else:
            i = randrange(N_SEG); j = randrange(N_SEG)
            var[0][i], var[1][j] = var[1][j], var[0][i]
            meta[0][i], meta[1][j] = meta[1][j], meta[0][i]
        counts[op] += 1

    write_genome(f"{OUTDIR}/genome_{label}.fa", var, label)
    # Record per-segment provenance + lengths after indels
    cumlen = [[0]*N_SEG for _ in range(N_CHR)]
    for c in range(N_CHR):
        running = 0
        for s in range(N_SEG):
            cumlen[c][s] = running
            running += meta[c][s]["len"]
    for c in range(N_CHR):
        for s in range(N_SEG):
            m = meta[c][s]
            provenance.append({
                "sample": label, "chr": f"chr{c+1}", "slot": s,
                "var_start": cumlen[c][s], "var_len": m["len"],
                "orig_chr": f"chr{m['orig_chr']+1}",
                "orig_seg": m["orig_seg"], "orient": m["orient"]
            })
    total_len = [sum(meta[c][s]["len"] for s in range(N_SEG)) for c in range(N_CHR)]
    print(f"{label}  subRate={div}  indelRate={indel_rate:.3f}  events={counts}  "
          f"chr_lengths={total_len}")

# Combined PanSN multi-fasta
samples = ["ref","div01","div02","div05","div10","div15"]
with open(f"{OUTDIR}/pangenome5.fa","w") as out:
    for s in samples:
        with open(f"{OUTDIR}/genome_{s}.fa") as fh:
            for L in fh:
                if L.startswith(">"):
                    m = re.match(r">(\S+)_(chr\d+)", L)
                    out.write(f">{m.group(1)}#1#{m.group(2)}\n")
                else:
                    out.write(L)

# Provenance TSV
with open(f"{OUTDIR}/truth.tsv","w") as fh:
    fh.write("sample\tchr\tslot\tvar_start\tvar_len\torig_chr\torig_seg\torient\n")
    # Reference rows
    for c in range(N_CHR):
        for s in range(N_SEG):
            fh.write(f"ref\tchr{c+1}\t{s}\t{s*SEG_LEN}\t{SEG_LEN}\tchr{c+1}\t{s}\t+\n")
    for r in provenance:
        fh.write("\t".join(str(r[k]) for k in
                 ("sample","chr","slot","var_start","var_len",
                  "orig_chr","orig_seg","orient")) + "\n")
print("wrote pangenome5.fa and truth.tsv")
