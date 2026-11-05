# The script behind each figure and table

Paths are relative to the repository root. Figure 2 of the main text is a schematic of the design and has no
generating script.

## Main figures

| Figure | Content | Script | Inputs |
|---|---|---|---|
| 1 | ADSP overview: participants, studies, sequencing centers, latent genetic classes | `07_figures_tables/figure1_overview.py` | panels from `07_figures_tables/figure1_overview_panels.R`; BIC curve from `01_preprocessing/fit_latent_genetic_classes.R` |
| 2 | The RIGA design | schematic | the cell grid of Supplementary Fig. 2 |
| 3 | Benchmark on observed ADSP genotypes | `07_figures_tables/figure3_benchmark.py` | panel tables from `07_figures_tables/figure3_benchmark_panel_tables.py`, which reads the outputs of `04_benchmark/single_variant_false_positive_counts.py`, `04_benchmark/gene_set_null_1000_scans_summary.py`, `04_benchmark/power_aggregate.py` and `04_benchmark/external_tools/` |
| 4 | Designs compared on observed data: single variants and gene sets | `07_figures_tables/figure4_figure5_design_comparison.py` | `07_figures_tables/figure4_extract_gene_set_results.R`, `figure4_extract_gene_comparison.R`, `figure4_extract_pooled_results.R` |
| 5 | Overlay of RIGA and pooled single-variant results | `07_figures_tables/figure4_figure5_design_comparison.py` | as Figure 4 |

## Supplementary figures

| Figure | Content | Script |
|---|---|---|
| 1 | Sequencing-center overview | `07_figures_tables/supplementary_figure1_sequencing_centers.R` |
| 2 | Project-by-center cell grid | `07_figures_tables/supplementary_figure2_cell_grid.py` |
| 3 | Rare-variant tools on the permuted null | `07_figures_tables/supplementary_figure3_rare_variant_tools.py` (inputs from `04_benchmark/gene_set_null_native_tests.py`, `rare_variant_tools_false_positives.py`) |
| 4 | Benchmark under ascertainment | `07_figures_tables/supplementary_figure4_benchmark_ascertainment.py --tag asc_clean_maf001` |
| 5 | Stratified designs against the reference | `07_figures_tables/supplementary_figure5_stratified_designs.py` (input from `04_benchmark/stratified_designs_compare.py`) |
| 6 | Synthetic multi-ancestry benchmark | `05_synthetic_benchmark/synthetic_benchmark_figure.py` on the runs of `synthetic_power_simulation.py` |
| 7 | Noncoding Manhattan and QQ plots | `07_figures_tables/supplementary_figure7_noncoding_manhattan_qq.py` (inputs from `03_scans/gene_set_primary_statistic.R`, `export_gene_centric_results.R`) |
| 8 | ncRNA comparison across designs | `07_figures_tables/figure4_figure5_design_comparison.py` |
| 9 | Burden power | `07_figures_tables/supplementary_figure9_burden_power.py` (uses `benchmark_figure_helpers.py`) |
| 10 | Partial-residual diagnostics of the null models | `07_figures_tables/supplementary_figure10_partial_residuals.py` from panels of `figure1_overview_panels.R` |

## Supplementary tables

The builders write CSV files; `07_figures_tables/supplementary_tables_specifications.py` holds the column
specifications and captions, `supplementary_tables_numbering.py` maps the published numbering to the analysis
numbering used inside the scripts, and `supplementary_tables_workbook.py` writes the workbook.

| Table | Content | Script |
|---|---|---|
| 1a to 1g | Data flow, self-reported groups, centers, projects, cohorts, the 79 cells, latent classes | `07_figures_tables/supplementary_table1_sample_description.R` |
| 2 | The 296 single variants at FDR 5% | `07_figures_tables/supplementary_table2_single_variants.R` |
| 3 | FinnGen R12 replication and GWAS Catalog audit | `06_replication/finngen_replication_and_catalog_audit.R` with `finngen_extract_summary_statistics.sh` |
| 4 | UK Biobank direct test of the non-APOE single variants | `06_replication/ukb/ukb_single_variant_direct_test.R`, `ukb_build_replication_tables.py` |
| 5a to 5c | Held-out halves | `02_design/split_sample_halves.R`, `split_half_null_models.R`, `03_scans/split_half_cross_replication.R` |
| 6, 6b, 6c | Calibration within self-reported groups | `02_design/ancestry_group_null_models.R`, `03_scans/calibration_within_ancestry_groups.R`, `consistency_across_ancestry_groups.R` |
| 7 | UK Biobank coding results for the named genes | `06_replication/ukb/ukb_coding_gene_lookup.R`, `ukb_snx22_recomputation.R` |
| 8 | Gene-by-mask sets at FDR 5% | `07_figures_tables/supplementary_table8_gene_sets.R` on `03_scans/gene_set_primary_statistic.R` |
| 9 | Carrier counts | `07_figures_tables/supplementary_table9_carrier_counts.R`, `supplementary_table9_carrier_counts_ncrna.R` |
| 10 | GWAS Catalog audit of the genes | `06_replication/finngen_replication_and_catalog_audit.R`, `gwas_catalog_search.py` |
| 11 | UK Biobank replication of the 14 gene sets | `06_replication/ukb/ukb_build_replication_tables.py` |
| 12 | Per-chromosome variant counts | `07_figures_tables/supplementary_table2_single_variants.R`, `supplementary_table12_variant_counts.R` |
| 13 | What each threshold admits | `07_figures_tables/supplementary_table8_gene_sets.R`, `supplementary_table13_threshold_tiers.R` |
| 14 | Software versions | `docs/ENVIRONMENT.md` |
| 15a | Single-variant false positives on the permuted null | `04_benchmark/single_variant_false_positive_intervals.py` |
| 15b | Gene-set null with and without the cMAC floor | `04_benchmark/gene_set_null_1000_scans_summary.py` |
| 16 | Diversity accounting | `03_scans/diversity_accounting.R` |
| 17 | Conditional analysis of the APOE region and TREM2 | `03_scans/conditional_analysis_apoe_trem2.R` |
| 18 | Robustness to the null-model specification | `03_scans/null_model_sensitivity_retest.R`, `null_model_sensitivity_summary.py` |
| 19 | Normal-approximation and saddlepoint tests of the 296 variants | `03_scans/single_variant_saddlepoint_retest.R`, `07_figures_tables/supplementary_table19_two_tests.py` |
| 20 | ADSP Gene Verification Committee catalog | compiled from topgenes.niagads.org |
| 21 | Scan comparisons across designs | `03_scans/gene_set_fdr_by_test.R`, `merge_single_variant_chunks.R` (counts per design) |
| 22 | Progressive-model ablation on the permuted null | `04_benchmark/benchmark_harness_model_ablation.py`, `model_ablation_summary.py` |
