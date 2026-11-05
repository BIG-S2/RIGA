# Software versions

All analyses ran on a Red Hat Enterprise Linux 9 cluster under SLURM (Supplementary Table 14 of the manuscript).

| Tool or package | Version | Role |
|---|---|---|
| R | 4.4.0 | analysis environment |
| STAAR | 0.9.8 | null model (`fit_nullmodel`), STAAR-O, saddlepoint burden test |
| STAARpipeline, STAARpipelineSummary | 0.9.8 | single-variant and gene-centric drivers, result summaries |
| SeqArray / SeqVarTools / gdsfmt / GDSArray | 1.46.3 / 1.44.0 / 1.42.1 / 1.26.0 | GDS input and output |
| GENESIS | 2.36.0 | mixed-model comparator |
| GMMAT | 1.4.2 | GLMM score test and SMMAT comparators |
| MetaSTAAR | 0.9.6.3 | exploratory ancestry-stratified meta-analyses only |
| mclust | 6.1.1 | Gaussian mixture (latent profile) analysis of PC1 to PC20 |
| MatchIt | 4.7.2 | matching for the sensitivity designs and the synthetic benchmark |
| Matrix, data.table | 1.7-4, 1.17.8 | sparse GRM representation, data handling |
| fastSparseGRM, KING | as distributed, 2.3.2 | unrelated set, PCs, sparse GRM |
| FAVORannotator | as distributed | functional annotation written into the annotated GDS |
| RUTH | 2021 release | ancestry-aware Hardy-Weinberg testing with 20 PCs |
| bcftools / htslib | 1.21 and 1.22 | genotype-level filtering, INFO recomputation, multi-allelic decomposition |
| PLINK, PLINK 2 | 1.90b7.4, 2.00a6LM | LD pruning (`--indep-pairwise 200 50 0.2`), format conversion |
| REGENIE | 3.2.9 | whole-genome regression comparator (Firth, leave-one-chromosome-out) |
| SAIGE-GENE+ | 1.3.3 | gene-based comparator, run from the authors' container image under Apptainer |
| Python (benchmark on observed genotypes) | 3.10; numpy 1.26.4, scipy 1.15.2, pandas 2.2.3 | `04_benchmark/`, figure scripts |
| Python (synthetic benchmark) | 3.12.7; numpy, scipy, pandas, scikit-learn, matplotlib | `05_synthetic_benchmark/` |

## Patched STAAR saddlepoint solver

The UK Biobank coding tables (Supplementary Tables 4, 7 and 11) were produced with STAAR and STAARpipeline 0.9.8
patched so that the Newton solver of the binary saddlepoint approximation verifies its root and falls back to a
bracketing bisection when it has not converged. `06_replication/staar_saddlepoint_patch/` holds the two patches, the
build script (`build_patched_staar.sh`, which fetches the upstream sources by commit, applies the patches and
installs them into a separate library), the unit test and validation scripts, and a README describing the defect and
the fix. Every other result in the paper used the unpatched 0.9.8 releases.
