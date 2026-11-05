#!/bin/bash
# v2 (2026-10-04): sparse GRM at relatedness cutoff 0.125 (env SAIGE_GRM). With the 0.05 GRM (19.4M entries) every pooled-design replicate whose
# variance component converged to a nonzero value (tau 0.03-0.08, ~35% of pooled replicates; 0% of stratified ones) never finished step 2.
# plasmode_rv_saige.sh -- SAIGE-GENE+ (official container, outputs/plasmode/rv_tools/containers/saige.sif) on one replicate of the
# plasmode permuted null: step 1 null GLMM with the SAIGE sparse GRM (createSparseGRM on the 281k panel; the pool is KING-unrelated
# so the GRM is near-diagonal), the covariate columns written by plasmode_rv_saige_inputs.py (full-rank numeric design columns, as
# GMMAT / GENESIS get) and a 3,000-marker step-1 plink file incl. MAC 1-20 variants for the categorical variance ratios
# (plasmode_rv_saige_step1_markers.sbatch); then step 2 group tests on the merged coding bed (plasmode_rv_merge.sh) with the
# SAIGE group file (plasmode_rv_groups.py).  2026-10-02 revision: annotation groups "nonsyn;syn,nonsyn,syn" (SAIGE separates
# annotations WITHIN a group by ';' and groups by ','; the earlier "syn,nonsyn;nonsyn;syn" produced the groups 'syn' and
# 'nonsyn;nonsyn;syn'); step-1 artifacts are kept in the log dir and reused; step 2 runs SAIGE_PAR (default 8) concurrent chunks
# of the group file (single-threaded per process; chr21 took 414 s, the genome would take ~11 h serially).
# Usage: plasmode_rv_saige.sh <replicate> <covset> <bfile prefix> <groups dir> <out root> [threads]
set -uo pipefail
REP=$1; COVSET=$2; BFILE=$3; GDIR=$4; OUTROOT=$5; THREADS=${6:-8}
PAR=${SAIGE_PAR:-8}
REPO=${ADSP_ROOT}
PY=${CONDA_ROOT}/envs/graph/bin/python
SIF=${REPO}/outputs/plasmode/rv_tools/containers/saige.sif
GRM=${SAIGE_GRM:-${REPO}/outputs/plasmode/rv_tools/saige_grm/pool_281k_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx}
STEP1_PLINK=${SAIGE_STEP1_PLINK:-${REPO}/outputs/plasmode/rv_tools/saige_step1/step1_markers}
SCR=${SIM_SCRATCH:-/dev/shm/plasmode_rv_saige}/${COVSET}_rep$(printf %03d "$REP"); mkdir -p "$SCR"
OUT=${OUTROOT}/saige; LOG=${OUT}/logs/${COVSET}_rep$(printf %03d "$REP"); mkdir -p "$LOG"
STEM=${OUT}/${COVSET}_rep$(printf %03d "$REP")
[ -s "${STEM}_genes.tsv.gz" ] && { echo "exists: ${STEM}_genes.tsv.gz"; exit 0; }
module load apptainer/1.5.0 2>/dev/null
RUN="apptainer exec -B /work,${HOME_DIR},/dev/shm ${SIF}"
echo "[saige] rep=${REP} covset=${COVSET} host=$(hostname) par=${PAR} $(date +%T)"
t0=$(date +%s)
if [ -s "${LOG}/step1.rda" ] && [ -s "${LOG}/step1.varianceRatio.txt" ]; then
  echo "[saige] step 1 reused from ${LOG}"; cp "${LOG}/step1.rda" "${LOG}/step1.varianceRatio.txt" "${SCR}/"
else
  $PY ${REPO}/revision/scripts/plasmode_rv_saige_inputs.py --replicate "$REP" --covset "$COVSET" --out "${SCR}/pheno.txt" || exit 2
  COVS=$(head -1 "${SCR}/pheno.txt" | tr '\t' '\n' | grep -vE "^(IID|Y)$" | paste -sd, -)
  $RUN step1_fitNULLGLMM.R --plinkFile=${STEP1_PLINK} --phenoFile="${SCR}/pheno.txt" --phenoCol=Y --covarColList="$COVS" \
    --sampleIDColinphenoFile=IID --traitType=binary --outputPrefix="${SCR}/step1" --nThreads="$THREADS" --LOCO=FALSE \
    --useSparseGRMtoFitNULL=TRUE --sparseGRMFile="$GRM" --sparseGRMSampleIDFile="${GRM}.sampleIDs.txt" \
    --isCateVarianceRatio=TRUE --useSparseGRMforVarRatio=TRUE --IsOverwriteVarianceRatioFile=TRUE > "${LOG}/step1.log" 2>&1 || { echo "step1 failed, see ${LOG}/step1.log"; tail -20 "${LOG}/step1.log"; exit 3; }
  cp "${SCR}/step1.rda" "${SCR}/step1.varianceRatio.txt" "${LOG}/"
fi
echo "[saige] step1 done in $(( $(date +%s) - t0 ))s"; t1=$(date +%s)
## step 2: PAR concurrent chunks of the group file (two lines per gene: var, anno; genes dealt round-robin)
$PY - "${GDIR}/saige.group" "${SCR}" "$PAR" <<'PYEOF'
import sys
src, scr, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(src) as f:
    lines = f.readlines()
assert len(lines) % 2 == 0 and all(l.split()[1] == "var" for l in lines[0::2]) and all(l.split()[1] == "anno" for l in lines[1::2])
fh = [open(f"{scr}/group_{i}.txt", "w") for i in range(k)]
for g in range(len(lines) // 2):
    fh[g % k].write(lines[2 * g]); fh[g % k].write(lines[2 * g + 1])
for h in fh:
    h.close()
print(f"[saige] {len(lines) // 2} genes -> {k} chunks", flush=True)
PYEOF
pids=(); nfail=0
for i in $(seq 0 $((PAR - 1))); do
  ( $RUN step2_SPAtests.R --bedFile=${BFILE}.bed --bimFile=${BFILE}.bim --famFile=${BFILE}.fam --AlleleOrder=alt-first \
      --SAIGEOutputFile="${SCR}/step2_${i}.txt" --LOCO=FALSE --GMMATmodelFile="${SCR}/step1.rda" --varianceRatioFile="${SCR}/step1.varianceRatio.txt" \
      --sparseGRMFile="$GRM" --sparseGRMSampleIDFile="${GRM}.sampleIDs.txt" --groupFile="${SCR}/group_${i}.txt" \
      --annotation_in_groupTest="nonsyn;syn,nonsyn,syn" --maxMAF_in_groupTest=0.01 --minMAF=0 --minMAC=0.5 \
      --is_output_markerList_in_groupTest=FALSE --is_single_in_groupTest=FALSE --MACCutoff_to_CollapseUltraRare=10 > "${LOG}/step2_${i}.log" 2>&1 ) &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || nfail=$((nfail + 1)); done
for i in $(seq 0 $((PAR - 1))); do [ -s "${SCR}/step2_${i}.txt" ] || { echo "chunk ${i}: no output"; nfail=$((nfail + 1)); }; done
[ "$nfail" -eq 0 ] || { echo "step2: ${nfail} chunk(s) failed, see ${LOG}/step2_*.log"; grep -il "error" ${LOG}/step2_*.log | head; exit 4; }
head -1 "${SCR}/step2_0.txt" > "${SCR}/step2.txt"
for i in $(seq 0 $((PAR - 1))); do tail -n +2 "${SCR}/step2_${i}.txt" >> "${SCR}/step2.txt"; done
echo "[saige] step2 done in $(( $(date +%s) - t1 ))s ($(( $(wc -l < "${SCR}/step2.txt") - 1 )) rows)"
$PY ${REPO}/revision/scripts/plasmode_rv_saige_parse.py --step2 "${SCR}/step2.txt" --out "${STEM}_genes.tsv.gz" --meta "${STEM}_meta.json" \
  --replicate "$REP" --covset "$COVSET" --step1-log "${LOG}/step1.log" --seconds-step1 $(( t1 - t0 )) --seconds-step2 $(( $(date +%s) - t1 )) || exit 5
head -3 "${SCR}/step2.txt" > "${LOG}/step2.head.txt"; cut -f2 "${SCR}/step2.txt" | sort | uniq -c > "${LOG}/step2.groups.txt"
rm -rf "$SCR"
echo "[saige] DONE $(date +%T)"
