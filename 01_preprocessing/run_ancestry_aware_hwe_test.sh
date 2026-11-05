#!/usr/bin/env bash
set -euo pipefail
chr=$1; start=$2; end=$3

VCF_DIR=${ADSP_ROOT}/processed
EIG=${ADSP_ROOT}/evec.intersect.vcforder.txt
SML=${ADSP_ROOT}/samples.intersect
OUT=${ADSP_ROOT}/data

vcf=${VCF_DIR}/r5_chr${chr}.split.vcf.gz
reg=chr${chr}:${start}-${end}

# # skip empty windows
# if ! bcftools view -H -r "$reg" "$vcf" | head -1 | grep -q .; then
#   echo "[skip] $vcf $reg – no variants" >&2
#   exit 0
# fi

tmp_bcf=$(mktemp --suffix=.bcf)
bcftools view -r "$reg" -S "$SML" -Ob -o "$tmp_bcf" "$vcf"
bcftools index -f "$tmp_bcf"

ruth --site-only \
     --vcf "$tmp_bcf" \
     --evec "$EIG" \
     --sm-list "$SML" \
     --field GT \
     --lambda 0 \
     --out "${OUT}/ruth.site_only.chr${chr}.${start}-${end}.vcf"

rm -f "$tmp_bcf" "$tmp_bcf".csi
