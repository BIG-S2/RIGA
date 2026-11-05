# Multi-Ancestry GWAS Simulation — Design Reference

This document records the cohort, ancestry, sequencing-centre, and phenotype
design used by the simulation pipeline in `synthetic_null_simulation.py` and
`synthetic_power_simulation.py`. The current settings are tuned to mimic the
correlation structure of real ADSP WGS data (51 cohorts × 15 sequencing
centres, n=38,757), as analysed in
[`src/ADSP_data_r5_seqcenter.Rmd`](../ADSP_data_r5_seqcenter.Rmd) and its
rendered output.

The design has evolved across versions; v7 is the current production
configuration.

---

## ADSP empirical targets (from ADSP_data_r5_seqcenter.html §4)

| Predictor → SeqCenter | Accuracy |
|---|---|
| PCs (PC1–20) only | 84.3% |
| Cohort only | 85.9% |
| PCs + Cohort | 98.2% |

Other features the simulation aims to reproduce:

- AD prevalence differs by ancestry (higher in AFR than in EUR cohorts).
- Sequencing centres process distinct cohort/ancestry mixes (Cohort×Centre
  table is 88.9% sparse, χ² ≪ 0.001).
- Each cohort has a primary centre but a long tail of secondary centres,
  so Cohort→Centre is high but not perfect.

---

## v7 (current production design)

### 1. Two simulated source populations

| | Pop A | Pop B |
|---|---|---|
| Reference data | population-specific simulated genotypes (LD-structured) | population-specific simulated genotypes (LD-structured) |
| Pool size | 10,000 | 10,000 |
| Baseline disease prevalence | **`POP_A_PREV = 0.15`** | **`POP_B_PREV = 0.08`** |

The per-population prevalence creates a generative Pop ↔ Y correlation
(higher AD risk in Pop A, mirroring AFR vs EUR ancestry-stratified AD
prevalence in ADSP).

### 2. Six studies (cohorts) with ancestry-skewed enrolment

| Study | Pop A draws | Pop B draws | role |
|:-:|:-:|:-:|---|
| 1 | 475 | 25  | Pop-A-skewed memory clinic |
| 2 | 475 | 25  | Pop-A-skewed biobank |
| 3 | 250 | 250 | multi-ancestry biobank |
| 4 | 250 | 250 | multi-ancestry clinic |
| 5 |  25 | 475 | Pop-B-skewed biobank |
| 6 |  25 | 475 | Pop-B-skewed clinic |

Constants in `synthetic_null_simulation.py`:

```python
N_COHORTS = 6
COHORT_POP_A_SIZES = [475, 475, 250, 250,  25,  25]
COHORT_POP_B_SIZES = [ 25,  25, 250, 250, 475, 475]
```

Total: 1500 Pop A + 1500 Pop B = **3000 subjects per replicate**. Study
predicts Pop with ~80% accuracy by design.

### 3. Per-pop, per-study case fractions (ascertainment + within-study Pop↔Y residual)

```python
COHORT_CASE_FRACTIONS_A = [0.85, 0.25, 0.55, 0.85, 0.25, 0.55]   # Pop A enriched
COHORT_CASE_FRACTIONS_B = [0.75, 0.15, 0.45, 0.75, 0.15, 0.45]   # Pop B (–0.10)
```

Two simultaneous confounds:

- **Across-study (ascertainment, Y↔Study)**: case fraction varies as
  `[case-enriched, control-enriched, balanced] × 2` regardless of Pop.
  Adjusting only for `study` neutralises this.
- **Within-study (ancestry, Y↔Pop | Study)**: case fraction is +0.05 for
  Pop A and –0.05 for Pop B in every study, so Y is *not* independent of
  Pop given study. Adjusting only for `study` leaves residual
  stratification — **this is what makes `combined_study` actually fail**
  (vs the v5 design where it was trivially calibrated).

Per-pop marginal case rates in the assembled cohort: **Pop A ≈ 0.595,
Pop B ≈ 0.405**.

### 4. Three sequencing centres with noisy Cohort→Centre map

```python
PRIMARY_STUDY_TO_CENTER = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3}
N_CENTERS = 3
CENTER_NOISE_PROB = 0.14
```

- Centre 1: studies 1, 2 → 95% Pop A → Pop-A-dominant centre
- Centre 2: studies 3, 4 → 50/50 → multi-ancestry centre
- Centre 3: studies 5, 6 → 95% Pop B → Pop-B-dominant centre

`_assign_centers` deterministically reassigns
`round(N × CENTER_NOISE_PROB) = 420` of 3000 subjects to a non-primary
centre. Empirical out-of-sample multinomial-logistic accuracies on a
single replicate match ADSP's targets:

| Predictor → Centre | v7 sim | ADSP |
|---|---|---|
| Study only | **0.852** | 0.859 |
| PCs (PC1–20) only | 0.41–0.67 (homogeneous / pure pops) <br> 0.64–0.93 (heterogeneous / admixed) | 0.843 |
| PCs + Study | 0.85–0.97 | 0.982 |

PCs→Centre under-shoots ADSP in pure-population/homogeneous-covariate
scenarios because the genotype simulation framework standardises within
population *before* concatenation, so PC1 is dominated by per-centre
batch shifts rather than allele-frequency divergence. Closing this gap
would require refactoring `genotype_simulation.py` to standardise across
the joint sample.

### 5. Sequencing-centre batch effects

Per-centre batch effects are added by `inject_batch_effects` after
cohort assembly:

- Each centre's batch-effect magnitude is sampled from
  `effect_size × U[0.5, 2.0]` (per-centre multiplier on the base
  `batch_size = 0.3`).
- Each centre's affected SNP fraction is sampled from
  `frac_affected × U[0.5, 1.5]` (base `batch_frac = 0.05` ⇒ ≈700 SNPs).
- Cross-centre heterogeneity is verified by `verify_batch_heterogeneity`
  (one-way ANOVA must show >10% of test SNPs significant).

### 6. Phenotype model

For null calibration runs (`synthetic_null_simulation.py`): logistic on
`age + gender` only, with `target_prev` set per-population.

For power runs (`synthetic_power_simulation.py`):
`logit P(Y=1) = β_age·age + β_gender·gender + Σ β_snp·X` where
`β_snp` is non-zero only at `num_causal` SNPs inside one randomly chosen
LD block.

Effect-size parameterisation:
```
β_causal = sqrt(h² · π² / (3 · K · (1 − h²)))
```
where `π²/3` is the logistic residual variance, giving an approximate
liability-scale heritability `h²`.

Genetic-architecture modes:

| Mode | Effects |
|---|---|
| `shared` | β_A = β_B (identical effects in both populations) |
| `pop_specific` | a fraction `pop_specific_frac` of causal variants are private to one pop (zero in the other) |
| `heterogeneous` | shared variants with bivariate-normal effects, cross-pop genetic correlation `rG` (default 0.8) |
| `mixed` | population-specific variants AND heterogeneous shared effects |

### 7. PCA & GMM ancestry assignment

- `sklearn PCA(svd_solver='randomized', n_components=20)` on the
  combined Xstd matrix.
- `assign_population_labels` fits a Gaussian mixture on
  `[PC1..PC20, one-hot(center)]` (centre dummies help the GMM
  disambiguate populations when ancestry alone is noisy).

### 8. Six analysis scenarios

Defined in `_scenario_configs`:

| Key | Pop A | Pop B | covariates |
|---|---|---|---|
| `2pure_homo` | AFR | EUR | homogeneous |
| `2pure_hetero` | AFR | EUR | heterogeneous |
| `1pure_1admixed_homo` | AFR | CEU5-YRI5 (admixed, multichr panel) | homogeneous |
| `1pure_1admixed_hetero` | AFR | CEU5-YRI5 | heterogeneous |
| `2admixed_homo` | CEU5-YRI5 | CEU5-CHS5 | homogeneous |
| `2admixed_hetero` | CEU5-YRI5 | CEU5-CHS5 | heterogeneous |

The `_multichr_simu_std.npy` admixed reference panels were rebuilt with
multi-chromosome LD references (chrs 1, 3–15, 18–22; 910K SNPs)
producing 14000-SNP LD-structured simulated genotypes, replacing the
original chr22-only versions which had insufficient ancestry-informative
markers to separate CEU5-YRI5 from CEU5-CHS5.

### 9. Six GWAS methods compared

All test `H₀: β_SNP = 0` via vectorised logistic score (Rao) test.

| Method | Covariates |
|---|---|
| `combined_study` | study + center + age + gender |
| `combined_pcs` | center + age + gender + PC1..PC20 |
| `combined_population` | assigned_pop + center + age + gender |
| `combined_study_pcs` | study + center + age + gender + PC1..PC20 |
| `combined_population_pcs` | assigned_pop + center + age + gender + PC1..PC20 |
| `matched_population` | greedy match on assigned_pop + center + gender + age (±tolerance); regression on age only |

Note: in v5 (3-study design) `study == center` exactly, so
`combined_pcs` and `combined_study_pcs` were the same regression. In v7
they are no longer collinear (`study` ≠ `center` due to the 6→3
collapse + noise).

### 10. Replicate count and significance threshold

The array is sized for `α = 5e-8` calibration:

```
M_null ≈ 12000 SNPs/run × 10000 runs = 1.2×10⁸ null tests
At α = 5e-8 → 6 expected FPs per method per scenario per config
At α = 5e-7 → 60 expected FPs (CV ≈ 13%)
At α = 5e-5 → 6000 expected FPs (CV ≈ 1%)
```

Single-replicate FPR estimates at α=5e-8 are pure noise; reliable
estimates require aggregation across all 10K replicates.

### 11. Output layout

```
outputs/power_runs_v{N}/
  power_run_{id:04d}_metrics.csv     # one row per scenario × config × method
  power_run_{id:04d}_pvals.npz       # raw M-dim p-value arrays
```

Aggregator: `src/simulation/aggregate_power_results.sh` (calls
`src/plot_power_figure.py`) → `outputs/power/`.

---

## Version history

| Version | Studies | Per-pop prev | Per-pop case frac | Centre map | Study↔Pop | Study→Ctr | Combined_study calibrated? |
|---|---|---|---|---|---|---|---|
| v3 | 6 mono-ancestry | 0.109 | n/a (1 pop per study) | identity | 1.0 (perfect) | n/a | complete-separation issues |
| v4 | 3 multi-ancestry, random | 0.109 | study-only | identity | ~0.5 | n/a | trivially calibrated |
| **v5** | 3 multi-ancestry, balanced | 0.109 | study-only (same per pop) | identity (study=center) | 0.5 | 1.000 | **calibrated** (no Y↔Pop \| Study) |
| v6 (cancelled) | 6 ancestry-skewed | 0.109 | study-only | 6→3, noise=0.002 | 0.80 | 0.998 | **still calibrated** (no within-study Y↔Pop) |
| **v7 (current)** | 6 ancestry-skewed | **0.15 / 0.08** | **per-pop differential (±0.05)** | 6→3, noise=0.14 | 0.80 | **0.852** | **fails honestly** ✓ |

### Key v7 design moves vs prior versions

1. **6 studies + 3 centres with non-trivial collapse.** Decouples
   `study` from `center` so PC-based methods (`combined_pcs`,
   `combined_study_pcs`, `combined_population_pcs`) become genuinely
   different regressions.
2. **Ancestry-skewed enrolment.** Studies 1–2 are 95% Pop A, 5–6 are 95%
   Pop B; introduces a real Cohort↔Ancestry correlation matching ADSP.
3. **Per-pop differential prevalence (0.15 vs 0.08).** Generative
   Pop↔Y correlation in the underlying logistic model.
4. **Per-pop differential case fractions within each study.**
   The crucial step: within every study, Pop A is +0.05 and Pop B is –0.05
   from the study's base case fraction. This breaks Y⊥Pop|Study and
   forces `combined_study` to leave residual stratification uncorrected.
5. **Centre noise calibrated to ADSP.** `CENTER_NOISE_PROB = 0.14` ⇒
   Study→Centre ≈ 0.852 (matching ADSP 0.859).

### Known limitations

- **PCs→Centre under-shoots ADSP** in pure-pop/homogeneous-covariate
  scenarios (~0.41–0.67 vs 0.84). Root cause: per-population genotype
  standardisation in `genotype_simulation.py` removes marginal
  allele-frequency differences across populations, so PC1 is driven by
  batch effects rather than ancestry. Heterogeneous-covariate and
  admixed scenarios reach 0.64–0.93 (the inhomogeneous covariate
  distribution restores PC structure).
- **Centre 2 is intentionally a "mixed" centre.** A more ADSP-like design
  with 15 centres of varied ancestry mixes would require expanding
  `PRIMARY_STUDY_TO_CENTER` and increasing `N_COHORTS` further.

---

## File map

| File | Role |
|---|---|
| `synthetic_null_simulation.py` | constants, cohort assembly, batch effects, PCA, GMM, GWAS methods |
| `synthetic_power_simulation.py` | causal-variant generation, per-pop preset-effect phenotypes, replicate runner |
| `run_power_array.sbatch` | SLURM array script (10000 reps × 7 scenarios × 27 configs × 6 methods) |
| `submit_power_array.sh` | wrapper that exports env vars and submits the array |
| `regenerate_admix_genotypes.sh` | rebuilds multichr admixed `.npy` files (256G mem, bigmem partition) |
| `aggregate_power_results.sh` | calls `plot_power_figure.py` to render `outputs/power/` summaries |
