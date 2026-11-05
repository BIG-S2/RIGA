# RIGA

Analysis code for **"Study design determines valid genetic discovery in heterogeneous sequencing consortia"**
(Xiaoqi Li, Xihao Li, Borui Tang, Tengfei Li, Weili Lin and Hongtu Zhu).

RIGA (Robust Integrative Genomic Association) analyzes a pooled sequencing consortium within the units through which
participants entered it. A study unit is one sequencing project processed at one sequencing center; every unit that
contains both cases and controls is retained, and participants are compared only within their unit. Applied to the
Alzheimer's Disease Sequencing Project (ADSP), this removes the false positives that pooled analysis produces when
the study of origin predicts both allele frequencies and case status, while keeping the power and the ancestral
diversity of the consortium. The repository holds the scripts behind every analysis, figure and table in the paper.
It contains no data.

## Contents

| Directory | Scripts |
|---|---|
| `01_preprocessing/` | `convert_vcf_to_gds.R`, `annotate_gds_with_favor.R`, `run_ancestry_aware_hwe_test.sh`, `variant_qc_config.R`, `build_variant_qc_label.R`, `extract_unrelated_set.R`, `compute_principal_components.R`, `merge_phenotype_and_pcs.R`, `fit_latent_genetic_classes.R` |
| `02_design/` | `prepare_sample_and_one_to_one_matching.R` (the sample preparation shared by all designs), `riga_study_units_null_model.R` (the RIGA design), the comparison designs `matched_design_null_model.R`, `pooled_design_null_model.R` and `pooled_with_project_term_null_model.R`, and `ancestry_group_null_models.R`, `split_sample_halves.R`, `split_half_null_models.R` |
| `03_scans/` | `single_variant_scan.R`, `gene_centric_coding_scan.R`, `gene_centric_noncoding_scan.R`, `gene_centric_ncrna_scan.R` (STAARpipeline drivers), `staar_o_ultra_rare_correction.R`, post-processing (`merge_single_variant_chunks.R`, `apply_variant_qc_to_single_variant_results.R`, `export_gene_centric_results.R`, `gene_set_primary_statistic.R`, `gene_set_fdr_by_test.R`, `merge_saddlepoint_rescan.R`, `single_variant_saddlepoint_retest.R`, `single_variant_qq_manhattan_by_design.R`), and the conditional, sensitivity, split-half and within-group analyses |
| `04_benchmark/` | the benchmark on observed ADSP genotypes: `benchmark_pool.py`, `permuted_phenotypes.py`, `benchmark_harness.py` and its common-variant, rare-burden, stratified-design and model-ablation variants, the burden tests, the false-positive counts and intervals, the 1,000-scan gene-set null, the rare-variant tools, and `external_tools/` (REGENIE, GMMAT, GENESIS on the same phenotypes) |
| `05_synthetic_benchmark/` | `synthetic_null_simulation.py`, `synthetic_power_simulation.py`, the comparator designs, aggregation and figure (`DESIGN.md` describes the generator) |
| `06_replication/` | FinnGen and GWAS Catalog scripts, `ukb/` (UK Biobank phenotype, null model, genome-wide coding scan, gene and variant tests, replication tables) and `staar_saddlepoint_patch/` (the patch to the STAAR saddlepoint solver used for the UK Biobank coding tables, with its build script, unit test and validation) |
| `07_figures_tables/` | one script per figure (`figure1_overview.py`, `figure3_benchmark.py`, `figure4_figure5_design_comparison.py`, `supplementary_figure1_...` to `supplementary_figure10_...`) and per supplementary table, plus the table specifications, numbering and workbook builders |
| `docs/` | `FIGURE_TABLE_MAP.md` (the script behind each figure and table), `GLOSSARY.md` (the working labels that appear inside the code), `ENVIRONMENT.md` (software versions) |

Figure 2 of the main text is a schematic of the design and has no script.

## Requirements

R 4.4.0 with STAAR 0.9.8, STAARpipeline 0.9.8, SeqArray, GENESIS 2.36.0, GMMAT 1.4.2, mclust 6.1.1, MatchIt 4.7.2
and data.table; Python 3.10 with numpy 1.26.4, scipy 1.15.2 and pandas 2.2.3 for the benchmark on observed
genotypes and Python 3.12 for the synthetic benchmark; bcftools 1.21/1.22, PLINK 1.9 and 2.0, RUTH, REGENIE 3.2.9
and SAIGE-GENE+ 1.3.3 as command-line tools. The complete list with versions is in `docs/ENVIRONMENT.md`.

## Data

ADSP Release 5 whole-genome sequencing and phenotype data (NIAGADS dataset NG00067.v16) are available to qualified
investigators through the NIAGADS Data Sharing Service after approval of a Data Use Certification
(https://dss.niagads.org/datasets/ng00067-v16/). UK Biobank data were used under application 22783. FinnGen R12
summary statistics are available from the FinnGen consortium. No individual-level data are redistributed here.

## Using the code
File directories were replaced by placeholders:
`${ADSP_ROOT}` is the analysis root that held the phenotype table (`pheno/`), the QC labels, PCs and sparse GRM
(`processed/`), the annotated GDS files, null models and scan outputs (`output/`) and the scripts themselves
(`src/`, `revision/scripts/`); `${WORK_DIR}`, `${HOME_DIR}`, `${CONDA_ROOT}` and `${R_LIBS_ROOT}` are the scratch
space, home directory, conda installation and R library root; `${FINNGEN_DIR}` holds the FinnGen summary
statistics. Replace them with your own locations, and point the few `source()` calls that name a script under
`${ADSP_ROOT}` at the copy in this repository. The directories are numbered in dependency order; the header
comment of each script states its inputs and outputs.

Two steps are not scripts in this repository. The LD pruning, KING kinship, unrelated set, PCA and sparse GRM were
produced with the fastSparseGRM package scripts (`getDivergence`, `extractUnrelated`, `runPCA` with 20 PCs and
`calcSparseGRM`) on the LD-pruned common-variant panel, and the STAARpipeline drivers were submitted per chromosome
with the standard STAARpipeline job templates.

## Citation and license

Please cite the manuscript (`CITATION.cff`). The code is released under the MIT License. The patches in
`06_replication/staar_saddlepoint_patch/` apply to STAAR and STAARpipeline, which are distributed by their authors
under their own terms.

## Funding and acknowledgments

This work was supported in part by the National Institute on Aging (NIA) of the National Institutes of Health (NIH)
under grants R01AG085581, RF1AG082938 and RF1AG098697 (to T.L. and H.Z.) and K01AG095286 (to T.L.), and by NIH grants
R01MH136055 (to H.Z. and T.L.), OT2OD038045 and R01AR082684 (to H.Z.), and R01DE024450 and R21HD120911 (to T.L.). The
funders had no role in study design, data collection and analysis, the decision to publish or the preparation of the
manuscript. This research has been conducted using the UK Biobank resource (application number 22783). We thank the
Alzheimer's Disease Sequencing Project and the NIAGADS Data Sharing Service for access to whole-genome sequencing
data, the participants and investigators of the contributing studies, and the University of North Carolina at Chapel
Hill Research Computing group for computational resources.

Contact: Hongtu Zhu (htzhu@email.unc.edu).
