#!/bin/bash
# Stage 2: pull the FinnGen R12 rows that fall in the ADSP hit windows.
# One streaming pass per endpoint, O(1) hash lookup per row on a 250 kb bin key.
set -uo pipefail

FG=${FINNGEN_DIR}/R12/summary_stats/release
WORK=${ADSP_ROOT}/revision/work/replication
BINS=$WORK/finngen_bins.txt      # 250 kb window bins, WINNER ARM ONLY
POSK=$WORK/finngen_pos.txt       # exact chr:pos of every hit, ALL ARMS
BIN=250000

# The AD / dementia endpoints, largest and most specific first.
# G6_ALZHEIMER is the primary: 13,964 cases vs 486,384 controls.
ENDPOINTS="G6_ALZHEIMER G6_AD_WIDE AD_LO_EXMORE F5_ALZHDEMENT KRA_PSY_DEMENTIA_EXMORE"

echo "window bins (winner arm): $(wc -l < "$BINS")   exact positions (all arms): $(wc -l < "$POSK")"

for EP in $ENDPOINTS; do
  f=$FG/finngen_R12_${EP}.gz
  o=$WORK/finngen_${EP}.tsv
  if [[ ! -f $f ]]; then echo "MISSING $f"; continue; fi
  echo -n "extracting $EP ... "
  zcat "$f" | awk -F'\t' -v OFS='\t' -v BIN=$BIN -v BINS="$BINS" -v POSK="$POSK" '
    BEGIN{ while((getline k < BINS) > 0) keepbin[k]=1
           while((getline k < POSK) > 0) keeppos[k]=1 }
    NR==1 { print "chrom","pos","ref","alt","rsids","nearest_genes",
                  "pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls",
                  "in_window"; next }
    { c=$1; sub(/^chr/,"",c);
      inw = ((c":"int($2/BIN)) in keepbin) ? 1 : 0
      # keep a row if it is inside a winner-arm window, OR sits within 5 bp of ANY
      # hit position in any arm (5 bp slack absorbs indel left-alignment differences)
      hit = 0
      for (d = -5; d <= 5 && !hit; d++) if ((c":"($2+d)) in keeppos) hit = 1
      if (inw || hit) print c,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,inw }
  ' > "$o"
  echo "$(( $(wc -l < "$o") - 1 )) rows"
done

# Also pull FinnGen's OWN genome-wide significant AD loci, for the reverse lookup
# (does ADSP recover known AD architecture?). Primary endpoint only.
echo -n "extracting G6_ALZHEIMER genome-wide hits ... "
zcat $FG/finngen_R12_G6_ALZHEIMER.gz | awk -F'\t' -v OFS='\t' '
  NR==1 { print "chrom","pos","ref","alt","rsids","nearest_genes",
                "pval","mlogp","beta","sebeta","af_alt"; next }
  $7+0 < 5e-8 { c=$1; sub(/^chr/,"",c); print c,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 }
' > $WORK/finngen_G6_ALZHEIMER_gwsig.tsv
echo "$(( $(wc -l < $WORK/finngen_G6_ALZHEIMER_gwsig.tsv) - 1 )) rows"

echo "DONE stage 2"
