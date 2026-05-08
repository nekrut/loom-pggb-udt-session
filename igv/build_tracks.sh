#!/usr/bin/env bash
set -euo pipefail
ENV=/media/anton/data/sandbox/udt1/.loom/env/bin
PATH=$ENV:$PATH
export PATH

cd "$(dirname "$0")"

REF=ref.fa
samtools faidx "$REF"

for V in div01 div02 div05 div10 div15; do
  IN=../genomes4/genome_${V}.fa
  # Rename variant headers from "div01_chr1" -> "${V}_chr1" already, fine.
  # Alignment BAM
  minimap2 -ax asm20 -t 8 "$REF" "$IN" 2>/dev/null \
    | samtools sort -@ 4 -o ${V}.bam -
  samtools index ${V}.bam

  # PAF + variants via paftools.js call
  minimap2 -cx asm20 --cs -t 8 "$REF" "$IN" 2>/dev/null \
    | sort -k6,6 -k8,8n \
    | k8 $ENV/paftools.js call -L 100 -l 100 -f "$REF" -s ${V} - \
    > ${V}.vcf
  bgzip -f ${V}.vcf
  tabix -p vcf ${V}.vcf.gz

  echo "[$V] BAM and VCF.gz built"
done
ls -la *.bam *.bai *.vcf.gz *.tbi
