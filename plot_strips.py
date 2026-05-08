#!/usr/bin/env python3
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns

truth = pd.read_csv("genomes3/truth.tsv", sep="\t")
graph = pd.read_csv("genomes3/graph_identity.tsv", sep="\t")

order = ["ref","div01","div02","div05","div10","div15"]

fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
for row, (df, label) in enumerate([(truth, "Truth (per-segment Hamming)"),
                                   (graph, "PGGB graph (node-sharing)")]):
    for col, chrom in enumerate(["chr1","chr2"]):
        ax = axes[row, col]
        sub = df[df["chr"]==chrom]
        sns.stripplot(data=sub, x="sample", y="identity",
                      order=order, ax=ax, jitter=0.25, size=7, alpha=0.85)
        sns.pointplot(data=sub, x="sample", y="identity", order=order,
                      ax=ax, color="black", linestyle="none", markers="_",
                      errorbar=None, markersize=18)
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"{label} — {chrom}")
        ax.set_xlabel("")
        ax.set_ylabel("identity to ref")
        ax.grid(axis="y", alpha=0.3)
plt.suptitle("Per-1kb-window identity vs. reference (10 windows × 6 samples)", y=1.0)
plt.tight_layout()
plt.savefig("genomes3/identity_strip.png", dpi=130, bbox_inches="tight")
print("wrote genomes3/identity_strip.png")
