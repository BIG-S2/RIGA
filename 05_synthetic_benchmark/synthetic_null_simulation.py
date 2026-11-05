"""
synthetic_null_simulation.py — Null simulation for multi-ancestry GWAS methods
============================================================================

Overview
--------
This module implements the **null simulation** (no causal variants) used to
evaluate type-I-error calibration of six GWAS association strategies across
seven ancestry-composition scenarios.  It is the companion to
``synthetic_power_simulation.py``, which adds causal variants to assess statistical power.

All exported symbols are imported by ``synthetic_power_simulation.py``; the two modules
share phenotype-simulation helpers, batch-effect injection, PCA, GMM
population labelling, and the six association-test wrappers.

Simulation pipeline (one replicate)
------------------------------------
For each of the six scenarios:

  1. **Genotype arrays** – Two pre-computed LD-structured standardised genotype
     matrices (M SNPs × N subjects) are loaded from disk.  M ≈ 14 000 SNPs
     spanning a single chromosome; N ≈ 10 000 subjects per ancestry group.
     See ``_load_default_arrays()`` for paths.

  2. **Phenotype simulation** – ``simulate_two_populations()`` draws age,
     gender, and binary disease outcome Y for each population under either a
     *homogeneous* covariate distribution (both populations identical) or a
     *heterogeneous* distribution (Pop A younger / more female, Pop B older /
     more male).  No genetic effect is included in the null simulation (h²=0).
     Disease prevalence is fixed at 10.9 %.

  3. **Cohort assembly** – ``sample_subpopulations()`` partitions each
     population into N_COHORTS=3 non-overlapping random groups and assigns
     study labels 1–3.  *Each study contains subjects from both Pop A and
     Pop B* so that study membership is not a proxy for ancestry (see the
     function docstring for the full design rationale and comparison with the
     earlier mono-ancestry design).

  4. **Batch effects** – ``inject_batch_effects()`` adds heterogeneous
     sequencing-centre artefacts to 5 % of SNPs, with a per-centre effect
     magnitude drawn from [0.5×, 2.0×] of a baseline shift of 0.3 standard
     deviations.  ``verify_batch_heterogeneity()`` asserts via one-way ANOVA
     that >10 % of tested SNPs show a significant centre difference.

  5. **PCA** – sklearn randomised PCA is fitted on the combined (M × N_total)
     genotype matrix transposed, yielding N_PCS=20 principal components per
     subject.  PCs capture both ancestry and batch variation.

  6. **GMM population labelling** – ``assign_population_labels()`` fits a
     Gaussian Mixture Model (n_components=2, full covariance, 20 initialisations)
     on the 20 PCs augmented with one-hot-encoded centre indicators.  The two
     GMM clusters are mapped to population labels by majority vote.  Accuracy
     relative to the true ``population`` column is printed and saved as
     ``pop_accuracy``.

  7. **Association testing** – ``run_all_methods_for_setting()`` runs all six
     methods (see below) on the full dataset, returning M-dimensional p-value
     arrays.

  8. **Metrics** – Lambda_GC (genomic inflation factor) and type-I error rates
     at α = 5×10⁻⁸ are computed per method and saved to a CSV file.

Six GWAS methods
-----------------
All methods test H0: β_SNP = 0 using a vectorised logistic score (Rao) test
implemented in ``assoc_Y_on_SNP_plus_covars()``.  They differ only in which
covariates are included in the null model:

  1. ``combined_study``
       Covariates: study + center + age + gender
       Controls for sequencing-centre batch effects and cohort membership, but
       NOT for fine-grained ancestry within each multi-ancestry cohort.
       Expected to show residual inflation when population stratification is
       strong.

  2. ``combined_pcs``
       Covariates: center + age + gender + PC1…PC20
       Uses the top 20 PCs to adjust for ancestry continuously.  Performs well
       when PCs cleanly capture population structure.

  3. ``combined_population``
       Covariates: assigned_pop (GMM label) + center + age + gender
       Uses discrete GMM-derived ancestry clusters.  Effective when GMM
       classification accuracy is high; degrades when populations overlap in
       PC space (see 2-admixed scenarios below).

  4. ``combined_study_pcs``
       Covariates: study + center + age + gender + PC1…PC20
       Combines cohort and PC adjustments.

  5. ``combined_population_pcs``
       Covariates: assigned_pop + center + age + gender + PC1…PC20
       Combines discrete population label and PC adjustments.

  6. ``matched_population``
       Case-control matching on assigned_pop + center + gender + age (±1 yr),
       then tests within matched pairs adjusting only for residual age
       difference.  Robust but discards unmatched subjects.

Six scenarios
---------------
  2pure_homo        AFR vs EUR, homogeneous age/gender
  2pure_hetero      AFR vs EUR, heterogeneous age/gender
  1pure_1admixed_homo   AFR vs CEU50%-YRI50%, homogeneous
  1pure_1admixed_hetero AFR vs CEU50%-YRI50%, heterogeneous
  2admixed_homo         CEU50%-YRI50% vs CEU50%-CHS50%, homogeneous
  2admixed_hetero       CEU50%-YRI50% vs CEU50%-CHS50%, heterogeneous

  Note on the 2-admixed scenarios (``2admixed_homo``, ``2admixed_hetero``)
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CEU5-YRI5 and CEU5-CHS5 both carry 50 % CEU ancestry and are simulated
  from chromosome 22 only.  Chr22 does not carry enough ancestry-informative
  markers to separate the non-CEU components (YRI vs CHS): the two admixed
  populations overlap extensively in PC1/PC2 space, GMM assignment accuracy
  approaches 50 % (chance), and all six methods therefore perform similarly.
  This is a genuine scientific finding — population structure that cannot be
  detected from the available genetic data cannot be corrected for — not a
  software defect.

Output files (per replicate)
-----------------------------
  run_{id:03d}_metrics.csv   One row per scenario × method, columns:
                               scenario_key, scenario_label, method,
                               n_valid_pvals, lambda_gc, type1_error_5e8
  run_{id:03d}_pvals.npz     Raw M-dimensional p-value arrays keyed by
                               "{scenario_key}__{method}".

Usage (command line)
--------------------
  python synthetic_null_simulation.py --replicate-id 1 --seed 42 \\
      --age-tolerance 1 --out-dir outputs/null_runs

SLURM array jobs: see run_simulation_array.sbatch.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import truncnorm, f_oneway, chi2
from scipy.linalg import qr as scipy_qr
import statsmodels.api as sm
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import argparse
import sys
from bisect import bisect_left
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.metrics import confusion_matrix, accuracy_score

_src_dir = Path(__file__).resolve().parent.parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))
from plot_simulation_main_figure import (
    compute_lambda_gc,
    type1_error_rate,
)

# Reproducibility/selection constants used across helpers.
RNG_SEED = 42
SUBSAMPLE_SIZE = 1000
N_PCS = 20

# Cohort design constants (used by sample_subpopulations).
#
# Real-ADSP-like design: 6 studies with **ancestry-skewed** enrolment.
# Each study draws from BOTH Pop A and Pop B but in unequal proportions,
# so `study` is partly predictive of ancestry (PCs ↔ study are correlated).
# Each study has its own case fraction (ascertainment confound, as before).
#
#                              PopA  PopB   case_frac   note
#  Study 1 (Pop-A clinic)      475   25     80%         AFR memory clinic
#  Study 2 (Pop-A biobank)     475   25     20%         AFR biobank
#  Study 3 (mixed centre 1)    250   250    50%         multi-ancestry biobank
#  Study 4 (mixed centre 2)    250   250    80%         multi-ancestry clinic
#  Study 5 (Pop-B biobank)     25    475    20%         EUR biobank
#  Study 6 (Pop-B clinic)      25    475    50%         EUR clinic (balanced)
#
# Total Pop A draws: 1500   Total Pop B draws: 1500   Total subjects: 3000
N_COHORTS = 6
COHORT_SIZE = SUBSAMPLE_SIZE // 2  # legacy; not used directly any more

# Per-study draws from each population (length must equal N_COHORTS).
COHORT_POP_A_SIZES = [475, 475, 250, 250, 25, 25]
COHORT_POP_B_SIZES = [25, 25, 250, 250, 475, 475]

# Differential case ascertainment per study (length must equal N_COHORTS).
# To create a Pop↔Y correlation that is NOT captured by `study`, the case
# fraction differs by ancestry *within* each study (the role: ancestry-
# stratified disease prevalence persists within every cohort, not just
# across cohorts).  Conditional on study, Y is therefore not independent
# of Pop, so adjusting only for `study` leaves residual stratification
# (which is what makes `combined_study` actually fail honestly).
#
# Pop A enrichment is +0.05 vs the base, Pop B is -0.05.
#
# SIM_INTERACTION_CONFOUNDING selects a three-valued interaction mode.  In
# EVERY mode each study's MEAN case fraction is identical, so the study main
# effect never moves and only the ANCESTRY effect is manipulated.
#
#   unset / "0" / "off"        no interaction; per-study A-B gap is +.10 x 6.
#   "1" / "center" / "centre"  ancestry x CENTRE interaction: the sign of the
#                              Pop A / Pop B enrichment is reversed in centre 2
#                              (studies 3 and 4).  Gaps +.10 +.10 -.10 -.10
#                              +.10 +.10; residual ancestry main effect +0.033.
#   "2" / "study"              ancestry x STUDY interaction: the sign is
#                              reversed in the SECOND study of every centre
#                              (studies 2, 4 and 6).  Gaps +.10 -.10 +.10 -.10
#                              +.10 -.10; residual ancestry main effect 0.000.
#
# The two interaction modes are ORTHOGONAL by construction, and that is the
# whole point of having both.  Decomposing the per-study gap into a
# between-centre and a within-centre component:
#
#   centre mode  var(gap) = 0.008889, 100% BETWEEN-centre -> matching on
#                (assigned_pop, center) cancels the entire interaction, because
#                both members of a matched pair carry the same centre term.
#   study  mode  var(gap) = 0.010000,   0% between-centre -> matching on
#                (assigned_pop, center) cancels NONE of it: within a centre the
#                two studies carry equal and opposite gaps that average to
#                exactly zero, so centre-level matching sees a null.  Only
#                matching on (assigned_pop, study) removes it.
#
# So the centre mode is favourable to RIGA-as-published by construction (the
# confounder was placed at RIGA's own matching resolution), and the study mode
# is the falsification test.  Reporting only the centre mode would be circular.
# `pvals_matched_study` exists so the study mode still has an arm that can
# succeed, which turns "RIGA fails" into the sharper and defensible statement
# that design-based control works exactly when the matching key contains the
# factor carrying the interaction.
#
# Total interaction magnitudes differ by 12.5% in variance between the two
# modes (0.008889 vs 0.010000).  If an exactly magnitude-matched study mode is
# ever wanted, use gaps of +/-0.09428 instead of +/-0.10.
#
# --- MECHANISM B: modes 3 and 4 -------------------------------------------
# Modes 1 and 2 act only through the case-fraction vectors, i.e. only through
# ANCESTRY.  Both were audited and both turned out to REDUCE the realized
# size-weighted marginal P(case|PopA) - P(case|PopB) rather than relocate it
# (0.1900 -> 0.1233 for the centre mode, -53% for the proposed study-level
# vectors), so "interaction present" was aliased with "less confounding to
# remove".  Worse, mode 2 has no genotype channel at all: `inject_batch_effects`
# shifts genotypes by CENTRE only, and the two studies inside every centre have
# identical ancestry composition, so E[G | pop, center, study] does not depend
# on study and a study-level ascertainment contrast cannot confound a test that
# already conditions on (assigned_pop, center).
#
# Modes 3 and 4 take the opposite route.  They leave the case-fraction vectors
# at their OFF-mode values -- so the realized ancestry marginal is baseline BY
# CONSTRUCTION, and the OFF mode is itself the matched-strength ancestry
# control -- and instead give `study` a genotype channel:
#
#   "3" / "studybatch"    a second, STUDY-resolution batch component is added
#                         on top of the unchanged centre component (see
#                         inject_group_batch_effects below).  Combined with the
#                         per-study case fractions [.80 .20 .50 .80 .20 .50],
#                         which are already part of the OFF design and are
#                         DISCORDANT within every centre, this makes
#                         E[G | Y, gender, assigned_pop, center] depend on Y.
#                         Invisible to a matching key or covariate set holding
#                         only `center`; removable by one holding `study`.
#   "4" / "pseudostudy"   the MATCHED-STRENGTH BATCH CONTROL for mode 3: the
#                         identical artefact (same affected-SNP counts, same
#                         effect sizes, same multiplier ladders, same group
#                         sizes, same nesting under centre) applied to a
#                         WITHIN-CENTRE PERMUTATION of the study labels.  The
#                         only thing destroyed is the alignment between the
#                         artefact and differential ascertainment.  Any
#                         difference between 3 and 4 is therefore attributable
#                         to that alignment, not to the amount of batch noise.
#
# WHAT MODES 3/4 DEMONSTRATE, stated plainly: this is BATCH x ASCERTAINMENT
# confounding, NOT ancestry confounding.  It is present at GMM assignment
# accuracy 1.0000, which is exactly why it escapes the tautology that matching
# on `assigned_pop` breaks any ancestry->Y path.  It answers "is design-based
# control robust to a confounder finer than its matching key", not "is
# design-based control better at removing ancestry confounding".
_INTERACT_MODES = {
    "": "off", "0": "off", "off": "off", "none": "off",
    "1": "center", "center": "center", "centre": "center",
    "2": "study", "study": "study",
    "3": "studybatch", "studybatch": "studybatch",
    "4": "pseudostudy", "pseudostudy": "pseudostudy",
    "5": "cohortmix", "cohortmix": "cohortmix",
    "6": "cohortmix_main", "cohortmix_main": "cohortmix_main",
    "7": "contanc", "contanc": "contanc",
    "8": "contanc_main", "contanc_main": "contanc_main",
}
# Modes that leave the phenotype side at its OFF-mode values and act only on
# the genotype side.  Kept as a set so the case-fraction branch below and the
# design-invariant assertions can treat them exactly like "off".
_GENO_CHANNEL_MODES = frozenset({"studybatch", "pseudostudy"})

# ---------------------------------------------------------------------------
# MECHANISM A -- modes 5 / 6 ("cohortmix", "cohortmix_main")
# ---------------------------------------------------------------------------
# THE DEFECT MECHANISM A REPAIRS.  In the OFF / centre / study designs the two
# studies inside every centre have IDENTICAL ancestry composition:
#
#     centre 1 = study1 (475,25) + study2 (475,25)
#     centre 2 = study3 (250,250) + study4 (250,250)
#     centre 3 = study5 (25,475) + study6 (25,475)
#
# so `study` carries NO information about ancestry beyond `center`, and a
# study-level ascertainment difference cannot confound a genotype test that
# already conditions on (assigned_pop, center).  That is why the study-level
# manipulation (mode 2) was falsified before it ran.
#
# MECHANISM A re-allocates the SAME 1500 Pop-A and 1500 Pop-B draws so that
# the two studies inside every centre have DIFFERENT ancestry mixes:
#
#     centre 1 = study1 (475,25)  + study2 (250,250)
#     centre 2 = study3 (250,250) + study4 (25,475)
#     centre 3 = study5 (25,475)  + study6 (475,25)
#
# The study -> centre map is UNCHANGED ({1:1,2:1,3:2,4:2,5:3,6:3}); only which
# study draws which ancestry mix moves.  Totals stay 1500 / 1500 / N=3000.
# This is the ADSP structure the real-data half of the revision established:
# USUHS holds 32 distinct cohorts with different ancestry composition, and the
# ancestry x Cohort interaction on case status survives conditioning on
# ancestry x SeqCenter.
#
#   MODE 5 "cohortmix"        the manipulation.  Per-study case fraction is
#                             IDENTICAL for Pop A and Pop B (gap d_s = 0 for
#                             EVERY study), so within a cohort ancestry is
#                             independent of Y.  The entire ancestry -> Y
#                             association is manufactured by cohort ancestry
#                             composition x cohort ascertainment:
#                               f_s = [.66 .50 .50 .34 .34 .66]
#                             A-rich studies (1, 6) are case-enriched, B-rich
#                             studies (4, 5) are control-enriched, balanced
#                             studies (2, 3) sit at 50%.
#
#   MODE 6 "cohortmix_main"   the MATCHED-STRENGTH MAIN-EFFECT CONTROL.  Same
#                             sizes, same centre-level ascertainment, same
#                             realized marginal, but the ancestry gap is a pure
#                             MAIN EFFECT: d_s = 0.16 for every study and the
#                             base rate is constant WITHIN each centre
#                             (f_c = [.55 .55 .45 .45 .50 .50]), so there is no
#                             ancestry x study and no ancestry x centre term.
#
# DESIGN ARITHMETIC (exact, on the literal vectors; the probe re-measures all
# of it on generated data).  Write the size-weighted marginal as
#     M = sum_s (nA_s - nB_s) f_s / 1500      [between-study composition term]
#       + sum_s (nA_s + nB_s) d_s / 3000      [within-study term]
# with f_s the per-study mean case fraction and d_s = fA_s - fB_s.
#
#                         M       R_centre   R_study   centre case fractions
#   off (baseline)      0.190      0.100      0.100     .545 / .650 / .305
#   cohortmix   (5)     0.192      0.166      0.000     .580 / .420 / .500
#   cohortmix_m (6)     0.190      0.160      0.160     .586 / .414 / .500
#
# where R_centre = sum_c W_c g_c / sum_c W_c with W_c = nA_c nB_c/(nA_c+nB_c)
# and g_c = P(case|A,c) - P(case|B,c)  (what a centre-conditional analysis
# sees), and R_study is the same statistic at study resolution with
# w_s = nA_s nB_s/(nA_s+nB_s)  (what a study-conditional analysis sees).
#
# So modes 5 and 6 are matched on the marginal (0.192 vs 0.190), on the
# centre-conditional residual (0.166 vs 0.160), and on the per-centre case
# fractions (.58/.42/.50 vs .586/.414/.50).  They differ in ONE thing: the
# study-conditional residual, 0.000 vs 0.160.  Mode 6 is therefore the control
# C1 demands -- any 5-vs-6 difference is attributable to WHERE the confounder
# sits, not to how much of it there is.  Because the two modes share the size
# vectors they also share GMM label informativeness, which the 5-vs-off
# comparison does NOT (see the probe's reported assignment accuracy).
#
# HONEST LIMIT, stated up front.  Genotype in this simulation depends only on
# (population panel, centre batch).  When `assigned_pop` equals true ancestry,
# matching on (assigned_pop, center) conditions on everything genotype depends
# on, and NO ancestry-side manipulation can touch it.  Mechanism A therefore
# can only bite through GMM label noise, i.e. in the admixed scenarios.  That
# is a real and reportable result ("RIGA is vulnerable exactly when its
# ancestry labels are noisy"), not a workaround, and the probe reports the
# admixed scenarios separately.
_COHORTMIX_MODES = frozenset({"cohortmix", "cohortmix_main"})

# ---------------------------------------------------------------------------
# MECHANISM C -- modes 7 / 8 ("contanc", "contanc_main")
# ---------------------------------------------------------------------------
# WHAT THE BRIEF ASKED FOR, AND WHY IT CANNOT BE BUILT LITERALLY.
# The brief asks for an interaction on a subject's CONTINUOUS ancestry
# coordinate (admixture proportion / position on the PC axis) rather than on
# the discrete `assigned_pop` label, so that two subjects sharing a label still
# differ in true ancestry and a matched pair does not cancel the confounder.
# Two facts about the generative model block that literally, and both were
# measured, not assumed (revision/logs/mechC_diag_66569791.out, and the
# per-panel statistics reprinted by revision/scripts/arenaA_mechC_verify.py):
#
#  F1  THERE IS NO PER-SUBJECT ANCESTRY COORDINATE.  Every panel is produced by
#      `simulate_genotypes(ld_sqrt_matrix=L, n_individuals=N)`, i.e. N i.i.d.
#      draws x_i = L z_i from ONE multivariate normal.  The "admixed" panels
#      (CEU5-YRI5, CEU5-CHS5) are not admixture simulations: they are ordinary
#      homogeneous panels whose LD reference happened to be an admixed sample.
#      No admixture proportion is drawn, retained, or retrievable.
#
#  F2  THE DISCRETE ANCESTRY LABEL HAS NO SINGLE-SNP GENOTYPE CHANNEL EITHER.
#      All four panels are standardised per SNP, so every panel has per-SNP
#      mean 0.000 and sd 1.000 to floating-point precision; the largest
#      between-panel per-SNP mean difference over 14,000 SNPs is 1.5e-16
#      (max |t| = 0.000 at n = 10,000/panel).  E[G_j | population] is therefore
#      IDENTICAL across populations for every SNP.  Populations differ only in
#      their COVARIANCE (which is what lets a full-covariance GMM separate them
#      at accuracy 1.0000 with no mean shift at all).  Consequently no
#      ancestry-side ascertainment manipulation -- discrete label or continuous
#      coordinate -- can bias a single-SNP score test, because the test reads
#      only E[G_j | Y].  This is a property of the benchmark, not of RIGA.
#
# WHAT MODES 7/8 DO INSTEAD, STATED WITHOUT EUPHEMISM.  They INSTALL the
# continuous coordinate that F1 says is missing, as a surrogate for fine-scale
# ancestry, and put the interaction on it:
#
#   a_i = u1' x_i,  u1 = leading eigenvector of the pooled panel covariance
#
# i.e. each subject's score on the first principal component of the genotype
# data -- the same quantity the real-data test used when it found PC1-5 x
# Cohort surviving conditioning on PC x SeqCenter (LRT p = 1.4e-37).  a_i has
# mean 0 in BOTH panels (F2), so it carries no between-population information
# whatever; it is purely WITHIN-label variation.  That is the point: it escapes
# the tautology that matching on `assigned_pop` breaks any ancestry->Y path,
# and it does so without relying on GMM label noise, so it is live in all six
# scenarios rather than only the two admixed ones.
#
# THE HONEST LABEL FOR THIS: it is a fine-scale WITHIN-POPULATION genotype
# gradient aligned with cohort recruitment.  It is what fine-scale ancestry
# looks like to a single-SNP test, and it is the closest thing to the brief's
# mechanism that this generative model can express.  It is NOT the simulation's
# own ancestry variable, because the simulation does not have one.
#
#   MODE 7 "contanc"       the manipulation.  Study s enrols subjects whose
#                          mean a is shifted by m_s = delta * [+1 -1 -1 +1 -1 +1],
#                          applied IDENTICALLY to that study's cases and its
#                          controls and identically to its Pop-A and Pop-B
#                          draws.  Within a study a is therefore independent of
#                          Y (nothing for study-level matching to miss), and
#                          sum_s m_s = 0 within EVERY centre (nothing for a
#                          between-centre correction to absorb).  The sign
#                          pattern is aligned with the OFF design's own per-
#                          study case fractions f_s = [.80 .20 .50 .80 .20 .50]:
#                          the case-enriched study of each centre pair gets +
#                          and the control-enriched one gets -.  Pooling within
#                          a centre then makes cases high-a and controls low-a
#                          -- a Simpson's-paradox confounder that is created
#                          entirely by cohort composition x cohort
#                          ascertainment.
#
#   MODE 8 "contanc_main"  the MATCHED-STRENGTH MAIN-EFFECT CONTROL.  The same
#                          coordinate, tilted by +delta' for cases and -delta'
#                          for controls in EVERY study identically -- a pure
#                          continuous-ancestry main effect with no study and no
#                          centre variation.  delta' is SOLVED (not chosen) so
#                          that the realized marginal E[a|case] - E[a|control]
#                          equals mode 7's.
#
# DESIGN ARITHMETIC (exact on the design vectors; the probe re-measures all of
# it on generated data).  Writing C_x = E[a|case, x] - E[a|control, x]:
#
#                   marginal C     C within centre     C within study
#   off                 0                0                   0
#   contanc      (7)   0.800 d         0.873 d              0
#   contanc_main (8)   0.800 d         0.800 d             0.800 d
#
# (measured on generated data at d = 0.5: 0.4005 / 0.4260 / 0.0020 for mode 7
# and 0.3956 / 0.3965 / 0.3964 for mode 8.)
#
# So 7 and 8 are matched exactly on the marginal, and mode 8 sits 8.3% BELOW
# mode 7 on the centre-conditional residual -- the quantity
# `matched_population` (RIGA as published) actually sees, so the control is if
# anything the WEAKER of the two by both measures.  They differ in exactly one
# number, the study-conditional residual: 0 against 0.800 d.  Any difference
# between them is therefore attributable to WHERE the confounder sits, never to
# how much of it there is.  Predicted: `matched_population` inflated equally in
# 7 and 8; `matched_study` calibrated in 7 and inflated in 8.
#
# The DISCRETE-label marginal P(case|PopA) - P(case|PopB) is untouched in both
# modes: they reuse the OFF case-fraction vectors verbatim, so C1 is satisfied
# at exactly the OFF baseline 0.1900 BY CONSTRUCTION, and OFF is itself the
# zero-strength control for the continuous channel.
#
# MEASURED OUTCOME (12 replicates x 6 scenarios x 3 modes on identical seeds;
# revision/tables/arenaA_mechC_*.csv, revision/logs/mechC_pilot_66579960_*.out).
# C1: the discrete-label marginal is +0.1900 with sd 0.0000 in ALL THREE modes,
# and the continuous-channel marginal is matched between 7 and 8 (0.4005 vs
# 0.3956).  The two modes differ in exactly one measured number, the
# within-study residual of a: +0.0020 +- 0.0042 in mode 7 against +0.3964 +-
# 0.0058 in mode 8.  C2, pooled, lambda_GC with its Monte-Carlo SE and FPR at
# 1e-3 as a multiple of nominal:
#
#                          off              contanc (7)        contanc_main (8)
#   matched_population   1.010  0.86x     1.252  30.7x         1.247  29.7x
#   matched_study        1.004  0.99x     0.988   0.94x        1.228  21.5x
#   combined_population  1.007  0.89x     1.289  40.0x         1.280  38.5x
#   combined_pcs         1.008  0.90x     1.064   3.3x         1.057   2.7x
#   combined_study_pcs   1.011  0.99x     1.016   1.05x        1.049   1.9x
#
# So RIGA as published is inflated 31x at alpha = 1e-3 while the same design
# matched one level finer is calibrated, and the mode-8 control shows this is
# not a magnitude effect: matched_population moves identically in 7 and 8
# (paired dlam +0.2420 +- 0.0134 vs +0.2366 +- 0.0112) while matched_study
# moves only in 8 (-0.0155 +- 0.0096 vs +0.2246 +- 0.0094).  It is also not an
# N effect: matched_study forms 705.8 pairs in mode 7 and 708.4 in mode 8
# against matched_population's 1026.6 / 1032.1 (0.69x in both) and is
# inflated 21x in one of them.  The excess false positives land where the
# mechanism put them -- in mode 7 matched_population's FPR is 280x nominal in
# the top decile of |u_j| and ~1.2x in the bottom seven, while matched_study is
# flat (0.68-1.14x) across every decile.
#
# The result holds in the four scenarios where GMM assignment accuracy is
# exactly 1.0000 (paired dlam for matched_population +0.2416 +- 0.0213) as
# strongly as in the two admixed ones (+0.2398 +- 0.0143).  That is the whole
# point of putting the interaction on a WITHIN-LABEL coordinate: unlike any
# manipulation acting through the ancestry LABEL, this one does not need label
# noise, so it is not confined to the admixed scenarios.
#
# delta IS NOT TUNED.  It is fixed at CONTANC_DELTA below, in units of the
# within-panel sd of a, and stated before the pilot ran.
_CONTANC_MODES = frozenset({"contanc", "contanc_main"})
_INTERACT_RAW = os.environ.get("SIM_INTERACTION_CONFOUNDING", "0").strip().lower()
if _INTERACT_RAW not in _INTERACT_MODES:
    raise ValueError(
        f"SIM_INTERACTION_CONFOUNDING={_INTERACT_RAW!r} is not one of "
        f"{sorted(set(_INTERACT_MODES))}"
    )
INTERACTION_MODE = _INTERACT_MODES[_INTERACT_RAW]
# Back-compat alias: True exactly when the old `== "1"` test was True.  Kept so
# that any external reader of `_INTERACT` keeps its historical meaning
# ("the centre-level interaction is on"), which is what the 1,500 finished
# replicates in outputs/power_runs_interact were generated under.
_INTERACT = INTERACTION_MODE == "center"

if INTERACTION_MODE == "center":
    COHORT_CASE_FRACTIONS_A = [0.85, 0.25, 0.45, 0.75, 0.25, 0.55]
    COHORT_CASE_FRACTIONS_B = [0.75, 0.15, 0.55, 0.85, 0.15, 0.45]
elif INTERACTION_MODE == "study":
    COHORT_CASE_FRACTIONS_A = [0.85, 0.15, 0.55, 0.75, 0.25, 0.45]
    COHORT_CASE_FRACTIONS_B = [0.75, 0.25, 0.45, 0.85, 0.15, 0.55]
elif INTERACTION_MODE == "cohortmix":
    # d_s = 0 for EVERY study: within a cohort, ancestry is independent of Y.
    COHORT_CASE_FRACTIONS_A = [0.66, 0.50, 0.50, 0.34, 0.34, 0.66]
    COHORT_CASE_FRACTIONS_B = [0.66, 0.50, 0.50, 0.34, 0.34, 0.66]
elif INTERACTION_MODE == "cohortmix_main":
    # d_s = 0.16 for every study, base rate constant within centre
    # (f_c = .55 .55 .45 .45 .50 .50): a pure ancestry main effect of the same
    # realized strength as mode 5.
    COHORT_CASE_FRACTIONS_A = [0.63, 0.63, 0.53, 0.53, 0.58, 0.58]
    COHORT_CASE_FRACTIONS_B = [0.47, 0.47, 0.37, 0.37, 0.42, 0.42]
else:
    COHORT_CASE_FRACTIONS_A = [0.85, 0.25, 0.55, 0.85, 0.25, 0.55]
    COHORT_CASE_FRACTIONS_B = [0.75, 0.15, 0.45, 0.75, 0.15, 0.45]

# Mechanism A re-allocates the ancestry mixes ACROSS studies so that the two
# studies inside every centre differ.  Totals are unchanged (1500 Pop A, 1500
# Pop B, N = 3000) and the study -> centre map is unchanged; only which study
# draws which mix moves.  Modes off / center / study / studybatch /
# pseudostudy never reach this branch, so their size vectors stay exactly as
# assigned above.
if INTERACTION_MODE in _COHORTMIX_MODES:
    COHORT_POP_A_SIZES = [475, 250, 250, 25, 25, 475]
    COHORT_POP_B_SIZES = [25, 250, 250, 475, 475, 25]
    assert sum(COHORT_POP_A_SIZES) == 1500 and sum(COHORT_POP_B_SIZES) == 1500, (
        f"mechanism-A sizes must still total 1500 / 1500, got "
        f"{sum(COHORT_POP_A_SIZES)} / {sum(COHORT_POP_B_SIZES)}"
    )
    # The whole point: the two studies inside a centre must have DIFFERENT
    # ancestry mixes, otherwise `study` has no channel `center` lacks.
    # _STUDY_TO_CENTER_DESIGN is defined a few lines below; a local copy is
    # used here and pinned against it by the assertion in that block.
    _S2C_MECHA = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3}
    for _c in sorted(set(_S2C_MECHA.values())):
        _members = [s for s, c in _S2C_MECHA.items() if c == _c]
        _mixes = {COHORT_POP_A_SIZES[s - 1] for s in _members}
        assert len(_mixes) == len(_members), (
            f"centre {_c}: studies {_members} have identical Pop-A draws "
            f"{_mixes}; mechanism A requires them to differ"
        )

# Aggregate per-study target (preserved for back-compat / logging only).
COHORT_CASE_FRACTIONS = [
    (a + b) / 2 for a, b in zip(COHORT_CASE_FRACTIONS_A, COHORT_CASE_FRACTIONS_B)
]

# --- Design invariants -------------------------------------------------------
# These run at import in every mode so the manipulation cannot silently drift.
# The study->centre map is spelled out here because PRIMARY_STUDY_TO_CENTER is
# defined much further down the file; a second assertion below its definition
# checks the two agree.
_STUDY_TO_CENTER_DESIGN = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3}
_GAP_TOL = 1e-12  # 0.85 - 0.75 is 0.09999999999999998 in binary FP, not 0.1.

_design_gaps = [a - b for a, b in zip(COHORT_CASE_FRACTIONS_A, COHORT_CASE_FRACTIONS_B)]
_design_centre_gap = {}
for _s, _c in _STUDY_TO_CENTER_DESIGN.items():
    _design_centre_gap.setdefault(_c, []).append(_design_gaps[_s - 1])
_design_centre_gap = {c: sum(v) / len(v) for c, v in _design_centre_gap.items()}
_design_net_gap = sum(_design_gaps) / len(_design_gaps)

# (i) STUDY MAIN EFFECT IS INVARIANT.  Holds in all three modes, and is exact in
#     binary floating point for these literals, so `==` is the right test.
assert len(COHORT_CASE_FRACTIONS_A) == len(COHORT_CASE_FRACTIONS_B) == N_COHORTS, (
    f"case-fraction vectors must have length N_COHORTS={N_COHORTS}"
)
#     Mechanism A (modes 5/6) is EXEMPT and says so out loud: it manufactures
#     the ancestry gap OUT OF the study main effect x cohort composition, so
#     the per-study case fractions necessarily move.  What is held instead is
#     the realized ancestry marginal and the per-centre case fractions -- both
#     asserted in the cohortmix branch below and re-measured on generated data
#     by revision/scripts/arenaA_cohortmix_probe.py.
assert INTERACTION_MODE in _COHORTMIX_MODES or (
    COHORT_CASE_FRACTIONS == [0.80, 0.20, 0.50, 0.80, 0.20, 0.50]), (
    f"mode {INTERACTION_MODE!r} changed the per-study MEAN case fractions to "
    f"{COHORT_CASE_FRACTIONS}; the interaction must move only the ancestry "
    f"effect, so any change here means the study main effect drifted too"
)
# (ii) and (iii) are mode-specific.
if INTERACTION_MODE == "study":
    # Each study's gap is +/-0.10, alternating in sign.
    assert all(abs(abs(g) - 0.10) < _GAP_TOL for g in _design_gaps), (
        f"study mode: every per-study A-B gap must be +/-0.10, got {_design_gaps}"
    )
    assert [g > 0 for g in _design_gaps] == [True, False, True, False, True, False], (
        f"study mode: gap signs must alternate within each centre, got {_design_gaps}"
    )
    # (ii) EVERY CENTRE'S MEAN ANCESTRY GAP IS ZERO.  This is what makes the
    #      manipulation invisible to centre-level matching.
    assert all(abs(v) < _GAP_TOL for v in _design_centre_gap.values()), (
        f"study mode: every centre's mean A-B gap must be 0, got "
        f"{_design_centre_gap}; if it is not, centre-level matching can absorb "
        f"part of the interaction and the experiment is no longer a clean test"
    )
    # (iii) NET ANCESTRY MAIN EFFECT IS ZERO.  This is what leaves nothing for
    #       an additive covariate or a PC to absorb.
    assert abs(_design_net_gap) < _GAP_TOL, (
        f"study mode: net ancestry main effect must be 0, got {_design_net_gap}; "
        f"a non-zero value would leave an additive component that a main-effect "
        f"correction could legitimately remove"
    )
elif INTERACTION_MODE == "center":
    assert all(abs(abs(g) - 0.10) < _GAP_TOL for g in _design_gaps), (
        f"centre mode: every per-study A-B gap must be +/-0.10, got {_design_gaps}"
    )
    # Sign reversal is confined to centre 2, so centre means are +.1/-.1/+.1 and
    # the net main effect is +1/3 of 0.10.
    assert abs(_design_centre_gap[1] - 0.10) < _GAP_TOL, _design_centre_gap
    assert abs(_design_centre_gap[2] + 0.10) < _GAP_TOL, _design_centre_gap
    assert abs(_design_centre_gap[3] - 0.10) < _GAP_TOL, _design_centre_gap
    assert abs(_design_net_gap - 0.10 / 3.0) < _GAP_TOL, (
        f"centre mode: net ancestry main effect must be +0.0333, got {_design_net_gap}"
    )
elif INTERACTION_MODE in _COHORTMIX_MODES:
    # ---- MECHANISM A design invariants -------------------------------------
    # Everything here is arithmetic on the literal design vectors.  It is a
    # guard against silent drift, NOT evidence: the probe re-derives all of it
    # from ACTUALLY GENERATED data, because "true by construction" is what the
    # last two designs also believed.
    _nA = COHORT_POP_A_SIZES
    _nB = COHORT_POP_B_SIZES
    _MECHA_TOL = 1e-9

    # (a) REALIZED SIZE-WEIGHTED MARGINAL P(case|PopA) - P(case|PopB) = 0.19
    #     +/- 0.01.  This is criterion C1, and it is the number the previous
    #     two designs got wrong by averaging the per-study gaps UNWEIGHTED.
    _mA = sum(n * f for n, f in zip(_nA, COHORT_CASE_FRACTIONS_A)) / sum(_nA)
    _mB = sum(n * f for n, f in zip(_nB, COHORT_CASE_FRACTIONS_B)) / sum(_nB)
    _mecha_marginal = _mA - _mB
    assert abs(_mecha_marginal - 0.19) <= 0.01, (
        f"mode {INTERACTION_MODE!r}: size-weighted marginal "
        f"P(case|A)-P(case|B) = {_mecha_marginal:.6f}, outside 0.19 +/- 0.01; "
        f"C1 (confounding held constant vs the OFF baseline) is violated"
    )

    # (b) WITHIN-STUDY ancestry-contrast-weighted residual, the quantity a
    #     study-adjusted / study-matched analysis actually sees.
    _w = [a * b / (a + b) for a, b in zip(_nA, _nB)]
    _mecha_R_study = sum(w * g for w, g in zip(_w, _design_gaps)) / sum(_w)

    # (c) CENTRE-conditional residual, the quantity a centre-adjusted /
    #     centre-matched analysis sees.
    _cen = sorted(set(_S2C_MECHA.values()))
    _nAc = {c: sum(_nA[s - 1] for s in _S2C_MECHA if _S2C_MECHA[s] == c) for c in _cen}
    _nBc = {c: sum(_nB[s - 1] for s in _S2C_MECHA if _S2C_MECHA[s] == c) for c in _cen}
    _gc = {}
    _fc = {}
    for c in _cen:
        _mem = [s for s in _S2C_MECHA if _S2C_MECHA[s] == c]
        _pA = sum(_nA[s - 1] * COHORT_CASE_FRACTIONS_A[s - 1] for s in _mem) / _nAc[c]
        _pB = sum(_nB[s - 1] * COHORT_CASE_FRACTIONS_B[s - 1] for s in _mem) / _nBc[c]
        _gc[c] = _pA - _pB
        _fc[c] = ((_nAc[c] * _pA + _nBc[c] * _pB) / (_nAc[c] + _nBc[c]))
    _W = {c: _nAc[c] * _nBc[c] / (_nAc[c] + _nBc[c]) for c in _cen}
    _mecha_R_centre = (sum(_W[c] * _gc[c] for c in _cen) / sum(_W.values()))

    if INTERACTION_MODE == "cohortmix":
        # The manipulation: ZERO within-study ancestry effect.  All of the
        # ancestry -> Y association is cohort composition x cohort
        # ascertainment, so a study-conditional analysis sees a null.
        assert all(abs(g) < _MECHA_TOL for g in _design_gaps), (
            f"cohortmix: every per-study A-B gap must be exactly 0, got "
            f"{_design_gaps}"
        )
        assert abs(_mecha_R_study) < _MECHA_TOL, _mecha_R_study
        # ... while a centre-conditional analysis does NOT.
        assert _mecha_R_centre > 0.10, (
            f"cohortmix: centre-conditional residual {_mecha_R_centre:.6f} is "
            f"no larger than the OFF baseline (0.10); the manipulation would "
            f"then have nothing for centre matching to miss"
        )
    else:  # cohortmix_main
        # The control: a PURE ancestry main effect of the same strength.  The
        # gap is identical in every study, hence in every centre.
        assert all(abs(g - 0.16) < _MECHA_TOL for g in _design_gaps), (
            f"cohortmix_main: every per-study A-B gap must be 0.16, got "
            f"{_design_gaps}"
        )
        assert all(abs(_gc[c] - 0.16) < _MECHA_TOL for c in _cen), (
            f"cohortmix_main: every centre's gap must be 0.16 (no ancestry x "
            f"centre term), got {_gc}"
        )
        assert abs(_mecha_R_study - 0.16) < _MECHA_TOL, _mecha_R_study

    print(
        f"[simulation] mechanism-A mode={INTERACTION_MODE!r} | sizes A={_nA} "
        f"B={_nB} | marginal={_mecha_marginal:.6f} | R_study="
        f"{_mecha_R_study:.6f} | R_centre={_mecha_R_centre:.6f} | per-centre "
        f"gaps { {c: round(v, 4) for c, v in _gc.items()} } | per-centre case "
        f"fractions { {c: round(v, 4) for c, v in _fc.items()} } | study "
        f"contrast weights { [round(x, 2) for x in _w] } "
        f"(share { [round(x / sum(_w), 3) for x in _w] })"
    )
else:
    # Covers "off", the Mechanism-B modes 3/4 AND the Mechanism-C modes 7/8,
    # all of which deliberately reuse the OFF case fractions so their realized
    # ancestry marginal equals baseline.
    assert (INTERACTION_MODE == "off"
            or INTERACTION_MODE in _GENO_CHANNEL_MODES
            or INTERACTION_MODE in _CONTANC_MODES), (
        f"unhandled INTERACTION_MODE {INTERACTION_MODE!r}"
    )
    assert all(abs(g - 0.10) < _GAP_TOL for g in _design_gaps), (
        f"off mode: every per-study A-B gap must be +0.10, got {_design_gaps}"
    )
    assert abs(_design_net_gap - 0.10) < _GAP_TOL, _design_net_gap

print(
    f"[simulation] SIM_INTERACTION_CONFOUNDING={_INTERACT_RAW!r} -> "
    f"mode={INTERACTION_MODE!r} | per-study gaps "
    f"{[round(g, 3) for g in _design_gaps]} | centre-mean gaps "
    f"{ {c: round(v, 3) for c, v in _design_centre_gap.items()} } | "
    f"net ancestry main effect {round(_design_net_gap, 4)}"
)

# Per-population baseline prevalence — creates a Pop↔Y correlation in the
# generative model.  Combined with the per-pop case fractions above,
# Pop A is genuinely higher-prevalence than Pop B, mirroring the ancestry-
# stratified late-onset AD prevalence in real ADSP.
POP_A_PREV = 0.15
POP_B_PREV = 0.08

# ---------------------------------------------------------------------------
# MECHANISM C machinery (modes 7 / 8).  Inert in every other mode: nothing
# below is called, no module constant above is changed, and no RNG is touched.
# ---------------------------------------------------------------------------
# delta, in units of the WITHIN-PANEL sd of the ancestry coordinate a.
# FIXED IN ADVANCE, NOT TUNED.  0.5 sd is the mean-a difference between the two
# cohorts inside a sequencing centre; real ADSP cohorts inside USUHS differ by
# far more than that on PC1, so this is a conservative setting.
#
# Two SPECIFICATION defects were found and fixed by the C1 probe before the
# pilot ran; neither was a change to delta and neither was made after seeing a
# calibration result.  (1) The coordinate was originally standardised POOLED
# across panels; measurement showed the leading pooled eigenvector carries ~100x
# more variance in one panel than the other, so a pooled-unit shift moved one
# population and not the other.  It is now standardised WITHIN each panel.
# (2) The tilt coefficient was originally used directly as the design
# parameter; because the case pool is ~10x more depleted than the control pool
# (prevalence 0.15 against per-study case fractions up to 0.85) a common
# coefficient saturated for cases and not for controls, leaking a within-study
# residual of -0.101 into the arm the design requires to be clean.  The design
# parameter is now the TARGET MEAN SHIFT and the coefficient is solved per pool.
CONTANC_DELTA = 0.5
# Sign of the per-study shift.  Aligned with the OFF design's own per-study
# case fractions f_s = [.80 .20 .50 .80 .20 .50] so that within each centre the
# case-enriched cohort is the high-a one:
#   centre 1 = (s1 f=.80 -> +) (s2 f=.20 -> -)
#   centre 2 = (s3 f=.50 -> -) (s4 f=.80 -> +)
#   centre 3 = (s5 f=.20 -> -) (s6 f=.50 -> +)
# Sums to zero over the two studies of EVERY centre, so there is no
# between-centre a gradient for a centre covariate or a PC to absorb.
CONTANC_STUDY_SIGN = [+1.0, -1.0, -1.0, +1.0, -1.0, +1.0]
# Deterministic power iteration for the ancestry axis: fixed seed and a fixed
# iteration count, so the axis is a property of the PANELS alone and does not
# vary by replicate or consume the caller's RNG.
CONTANC_AXIS_SEED = 20260902
CONTANC_AXIS_ITERS = 60
_CONTANC_AXIS_CACHE: Dict[tuple, tuple] = {}


def _contanc_expected_C(m_s: List[float]) -> float:
    """Expected E[a|case] - E[a|control], pooled over both populations.

    Exact given the design vectors, for a per-study mean shift ``m_s`` applied
    identically to that study's cases and controls.  This is the realized
    strength the mode-8 control has to reproduce; it is SOLVED for, not chosen.
    """
    num_c = num_k = den_c = den_k = 0.0
    for s in range(N_COHORTS):
        for n, f in ((COHORT_POP_A_SIZES[s], COHORT_CASE_FRACTIONS_A[s]),
                     (COHORT_POP_B_SIZES[s], COHORT_CASE_FRACTIONS_B[s])):
            num_c += n * f * m_s[s]
            den_c += n * f
            num_k += n * (1.0 - f) * m_s[s]
            den_k += n * (1.0 - f)
    return num_c / den_c - num_k / den_k


def _contanc_residual(m_s: List[float], group_of_study: Dict[int, int]) -> float:
    """Case/control-contrast-weighted residual of a within a grouping.

    The continuous-coordinate analogue of the ancestry statistic
    ``sum_s w_s g_s / sum_s w_s``: groups are weighted by
    ``v_g = n_case,g * n_ctrl,g / n_g`` (the within-group case/control contrast)
    and ``C_g = E[a|case,g] - E[a|control,g]``.  Pass the identity map for study
    resolution and PRIMARY_STUDY_TO_CENTER for centre resolution.
    """
    acc: Dict[int, List[float]] = {}
    for s in range(N_COHORTS):
        g = group_of_study[s + 1]
        d = acc.setdefault(g, [0.0, 0.0, 0.0, 0.0])   # nc, nk, sum a*nc, sum a*nk
        for n, f in ((COHORT_POP_A_SIZES[s], COHORT_CASE_FRACTIONS_A[s]),
                     (COHORT_POP_B_SIZES[s], COHORT_CASE_FRACTIONS_B[s])):
            d[0] += n * f
            d[1] += n * (1.0 - f)
            d[2] += n * f * m_s[s]
            d[3] += n * (1.0 - f) * m_s[s]
    num = den = 0.0
    for nc, nk, sc, sk in acc.values():
        if nc <= 0 or nk <= 0:
            continue
        v = nc * nk / (nc + nk)
        num += v * (sc / nc - sk / nk)
        den += v
    return num / den if den > 0 else 0.0


# The mode-7 design vector and the mode-8 tilt that matches its realized
# strength.  Computed at import in every mode (they are pure arithmetic on the
# design constants and change nothing), so the probe can read them without
# re-deriving them.
CONTANC_M_S = [CONTANC_DELTA * s for s in CONTANC_STUDY_SIGN]
CONTANC_C7 = _contanc_expected_C(CONTANC_M_S)
# Mode 8 shifts cases by +d' and controls by -d' everywhere, so its marginal C
# is exactly 2*d'.  Solve, do not choose.
CONTANC_MAIN_DELTA = CONTANC_C7 / 2.0
_CONTANC_IDENT = {s: s for s in range(1, N_COHORTS + 1)}

if INTERACTION_MODE in _CONTANC_MODES:
    # Design invariants for Mechanism C, asserted at import so the manipulation
    # cannot drift silently.  As with Mechanism A these are guards, not
    # evidence: revision/scripts/arenaA_mechC_verify.py re-measures every one of
    # them on ACTUALLY GENERATED data.
    _CONTANC_TOL = 1e-9
    # (a) The DISCRETE-label channel is untouched -- the case fractions are the
    #     OFF vectors verbatim, so the realized size-weighted marginal
    #     P(case|PopA) - P(case|PopB) is the 0.1900 baseline exactly.  C1 holds
    #     by construction and OFF is the zero-strength control.
    assert COHORT_CASE_FRACTIONS_A == [0.85, 0.25, 0.55, 0.85, 0.25, 0.55], \
        COHORT_CASE_FRACTIONS_A
    assert COHORT_CASE_FRACTIONS_B == [0.75, 0.15, 0.45, 0.75, 0.15, 0.45], \
        COHORT_CASE_FRACTIONS_B
    assert COHORT_POP_A_SIZES == [475, 475, 250, 250, 25, 25], COHORT_POP_A_SIZES
    assert COHORT_POP_B_SIZES == [25, 25, 250, 250, 475, 475], COHORT_POP_B_SIZES
    # (b) No between-centre a gradient: the two studies of every centre carry
    #     equal and opposite shifts and equal sizes.
    for _c in sorted(set(_STUDY_TO_CENTER_DESIGN.values())):
        _mem = [s for s, c in _STUDY_TO_CENTER_DESIGN.items() if c == _c]
        assert abs(sum(CONTANC_M_S[s - 1] for s in _mem)) < _CONTANC_TOL, (
            f"centre {_c}: per-study a-shifts {[CONTANC_M_S[s - 1] for s in _mem]} "
            f"do not cancel; a between-centre gradient would be absorbable by "
            f"the centre covariate and by the PCs"
        )
        assert len({COHORT_POP_A_SIZES[s - 1] + COHORT_POP_B_SIZES[s - 1]
                    for s in _mem}) == 1, "unequal study sizes break the cancellation"
    # (c) The sign is aligned with the per-study case fraction inside each
    #     centre -- this alignment IS the confounder.
    for _c in sorted(set(_STUDY_TO_CENTER_DESIGN.values())):
        _mem = sorted(s for s, c in _STUDY_TO_CENTER_DESIGN.items() if c == _c)
        _hi = max(_mem, key=lambda s: COHORT_CASE_FRACTIONS[s - 1])
        assert CONTANC_M_S[_hi - 1] > 0, (
            f"centre {_c}: the case-enriched study {_hi} must carry the +delta "
            f"shift, else the manipulation cancels instead of confounding"
        )
    # (d) 7 and 8 are matched on the marginal and on the centre-conditional
    #     residual, and differ only on the study-conditional residual.
    # `_STUDY_TO_CENTER_DESIGN` is the local copy of the study -> centre map
    # (PRIMARY_STUDY_TO_CENTER is defined further down and is pinned equal to
    # it by an assertion at that point).
    _R7_study = _contanc_residual(CONTANC_M_S, _CONTANC_IDENT)
    _R7_centre = _contanc_residual(CONTANC_M_S, _STUDY_TO_CENTER_DESIGN)
    _R8_study = _contanc_residual([0.0] * N_COHORTS, _CONTANC_IDENT) + 2 * CONTANC_MAIN_DELTA
    assert abs(_R7_study) < _CONTANC_TOL, (
        f"contanc: within-study residual of a must be exactly 0, got "
        f"{_R7_study}; a non-zero value would give study-level matching "
        f"something to remove and destroy the contrast with mode 8"
    )
    # Modes 7 and 8 must be MATCHED on the confounding magnitude, or a 7-vs-8
    # difference is a magnitude difference and the experiment proves nothing --
    # exactly the defect that killed the previous two designs.  A one-parameter
    # control cannot match BOTH the marginal and the centre-conditional
    # residual exactly, so delta' is solved on the MARGINAL (mode 8's tilt is
    # +/-delta' everywhere, hence its marginal C is exactly 2*delta' = C7) and
    # the centre-conditional residual then lands ~8% BELOW mode 7's.  That is
    # the conservative direction: the control has less-or-equal confounding by
    # BOTH measures, so if it nonetheless inflates an arm that mode 7 leaves
    # calibrated, the difference cannot be blamed on magnitude.  10% is the
    # asserted tolerance and both numbers are printed.
    _R8_centre = 2 * CONTANC_MAIN_DELTA
    assert _R8_centre <= _R7_centre + 1e-12, (
        f"contanc: mode-8 centre-conditional residual {_R8_centre:.6f} EXCEEDS "
        f"mode 7's {_R7_centre:.6f}; the control must not carry more "
        f"confounding than the manipulation it controls for"
    )
    assert abs(_R7_centre - _R8_centre) <= 0.10 * abs(_R7_centre), (
        f"contanc: centre-conditional residual of a is {_R7_centre:.6f} in "
        f"mode 7 but {_R8_centre:.6f} in mode 8, more than 10% apart; the two "
        f"modes are no longer matched on what `matched_population` sees"
    )
    print(
        f"[simulation] mechanism-C mode={INTERACTION_MODE!r} | delta="
        f"{CONTANC_DELTA} sd | per-study a-shift m_s={CONTANC_M_S} | "
        f"expected marginal C=E[a|case]-E[a|ctrl]={CONTANC_C7:.6f} | "
        f"mode-8 tilt delta'={CONTANC_MAIN_DELTA:.6f} | residual of a: "
        f"study {_R7_study:.3e} (mode 8: {_R8_study:.6f}), centre "
        f"{_R7_centre:.6f} (mode 8: {_R8_centre:.6f}, "
        f"{100 * (_R8_centre / _R7_centre - 1):+.1f}%)"
    )


def contanc_ancestry_axis(X_A: np.ndarray, X_B: np.ndarray) -> np.ndarray:
    """Leading eigenvector of the POOLED panel covariance, by power iteration.

    Returns a unit vector ``u`` of length M.  Deterministic (fixed seed, fixed
    iteration count) and cached per panel pair, so the axis is a property of
    the panels alone: it does not vary across replicates and never touches the
    caller's RNG.

    Power iteration is used rather than a full SVD because only the leading
    direction is needed and ``X_A X_A' + X_B X_B'`` is 14,000 x 14,000; four
    mat-vecs per iteration cost ~1e8 flops and the 14,000-square matrix is
    never formed.
    """
    key = (id(X_A), X_A.shape, id(X_B), X_B.shape)
    hit = _CONTANC_AXIS_CACHE.get(key)
    if hit is not None:
        return hit[0]
    M = X_A.shape[0]
    if X_B.shape[0] != M:
        raise ValueError(f"panels have different SNP counts: {X_A.shape} vs {X_B.shape}")
    rs = np.random.default_rng(CONTANC_AXIS_SEED)
    u = rs.normal(size=M)
    u /= np.linalg.norm(u)
    lam = 0.0
    for _ in range(CONTANC_AXIS_ITERS):
        v = X_A @ (X_A.T @ u) + X_B @ (X_B.T @ u)
        lam = float(np.linalg.norm(v))
        if lam <= 0:
            raise RuntimeError("degenerate ancestry axis (zero eigenvalue)")
        u = v / lam
    # Sign is arbitrary for an eigenvector; pin it so runs are comparable.
    if u[np.argmax(np.abs(u))] < 0:
        u = -u
    _CONTANC_AXIS_CACHE[key] = (u, X_A, X_B)      # hold refs so the ids stay live
    return u


def contanc_scores(X_A: np.ndarray, X_B: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Per-subject continuous ancestry coordinate for both full panels.

    ``a_i = u' x_i``, standardised WITHIN EACH PANEL, so CONTANC_DELTA is read
    as "half a within-population sd of the ancestry coordinate" and means the
    same thing in both populations.

    Pooled standardisation was tried first and rejected on measurement, not on
    taste.  The panels have identical per-SNP means but different covariances,
    and the leading pooled eigenvector turns out to be a direction in which one
    panel carries ~100x the variance of the other (measured: sd 0.143 vs 1.407
    in pooled units for AFR/EUR).  A common tilt in pooled units therefore
    moves one population and leaves the other essentially untouched.  Per-panel
    standardisation makes the RELATIVE displacement equal; the ABSOLUTE
    genotype consequence still differs by the ratio of the raw sds, and
    ``contanc_scores_info`` reports that ratio so it is never hidden.

    Both panels have per-SNP mean 0 by construction, so ``a`` has mean 0 in
    EACH panel either way: the coordinate carries no between-population
    information and every shift imposed on it is WITHIN-label variation.
    """
    aA, aB, _ = contanc_scores_info(X_A, X_B)
    return aA, aB


def contanc_scores_info(X_A: np.ndarray, X_B: np.ndarray):
    """``contanc_scores`` plus the raw scale of each panel's coordinate.

    Returns ``(aA_z, aB_z, info)`` where ``info`` carries the RAW (pre-
    standardisation) sd of a in each panel.  A CONTANC_DELTA shift of the
    standardised coordinate moves SNP j by ``delta * u_j * sd_raw`` in sd
    units, so these two numbers are what decide how hard the mechanism hits
    each population.
    """
    u = contanc_ancestry_axis(X_A, X_B)
    aA = np.asarray(X_A.T @ u, dtype=np.float64)
    aB = np.asarray(X_B.T @ u, dtype=np.float64)
    sA, sB = float(aA.std()), float(aB.std())
    if sA <= 0 or sB <= 0:
        raise RuntimeError("degenerate ancestry coordinate (zero sd in a panel)")
    info = {"sd_raw_A": sA, "sd_raw_B": sB, "sd_ratio": max(sA, sB) / min(sA, sB)}
    return (aA - aA.mean()) / sA, (aB - aB.mean()) / sB, info


# --- calibrating the enrolment tilt ----------------------------------------
# The tilt is exponential: a study enrols from its pool with probability
# proportional to exp(c * a), drawn without replacement by the Gumbel-top-k
# identity.  For an infinite pool and a standard-normal coordinate that shifts
# the selected mean by exactly c, but the pools here are NOT effectively
# infinite: with prevalence 0.15 the design needs ~893 of the ~1500 available
# Pop-A cases and only ~608 of the ~8500 available Pop-A controls, so at a
# common c the case draw saturates and the control draw does not.  Running the
# uncalibrated version measured exactly that: the within-study residual of a,
# which the design requires to be 0, came out at -0.101 against a centre-
# conditional residual of +0.289 -- a 35% leak straight into the arm the
# experiment is trying to keep clean.
#
# So c is not used as the design parameter.  The design parameter is the TARGET
# MEAN SHIFT m, and c is solved for, separately for every (study, population,
# case/control) pool, against that pool's actual size and composition.  The
# solve uses the standard conditional-Poisson approximation to successive
# sampling: inclusion probability pi_i = 1 - exp(-t * exp(c*a_i)) with t fixed
# by sum_i pi_i = k, giving E[mean of selected a] = sum_i pi_i a_i / k.  Both
# solves are monotone, so plain bisection is enough and the result is
# deterministic (no RNG is consumed).
CONTANC_SOLVE_ITERS = 40
CONTANC_C_BOUND = 12.0


def _contanc_expected_mean(a: np.ndarray, k: int, c: float) -> float:
    """E[mean of a over a Gumbel-top-k draw of size k at tilt c]."""
    # Normalise the weights to max 1.  Subtracting `a.max()` would only do
    # that for c > 0; half the studies carry a NEGATIVE target and hence a
    # negative c, for which `exp(c*(a - a.max()))` grows instead of shrinking
    # and can leave the required t outside any fixed bisection window.
    ca = c * a
    w = np.exp(ca - ca.max())
    lo, hi = 1e-30, 1e30
    for _ in range(CONTANC_SOLVE_ITERS):
        t = np.sqrt(lo * hi)
        s = float(np.sum(-np.expm1(-t * w)))
        if s < k:
            lo = t
        else:
            hi = t
    t = np.sqrt(lo * hi)
    pi = -np.expm1(-t * w)
    tot = float(pi.sum())
    return float(np.dot(pi, a) / tot) if tot > 0 else float(a.mean())


def contanc_tilt_coef(a: np.ndarray, k: int, target: float) -> Tuple[float, float]:
    """Solve for the tilt ``c`` whose expected selected mean of ``a`` is ``target``.

    Returns ``(c, achieved)``.  ``achieved`` is the expected mean at the solved
    ``c``; when the target is unreachable because the pool is too depleted the
    bisection lands on the bound and ``achieved`` reports how far it actually
    got, so saturation is visible rather than silent.
    """
    if k <= 0 or k >= a.size or target == 0.0:
        return 0.0, float(a.mean()) if k > 0 else 0.0
    lo, hi = -CONTANC_C_BOUND, CONTANC_C_BOUND
    for _ in range(CONTANC_SOLVE_ITERS):
        c = 0.5 * (lo + hi)
        if _contanc_expected_mean(a, k, c) < target:
            lo = c
        else:
            hi = c
    c = 0.5 * (lo + hi)
    return c, _contanc_expected_mean(a, k, c)


def contanc_study_tilts(study_1based: int) -> Tuple[float, float]:
    """(m_case, m_ctrl) TARGET mean shifts of a for one study, in sd units.

    Mode 7 shifts a study's cases and controls IDENTICALLY (so a is independent
    of Y within the study); mode 8 shifts cases up and controls down by the
    same amount in EVERY study (a pure main effect).  These are target MEANS,
    not tilt coefficients: the coefficient that delivers each one is solved for
    against the actual pool by ``contanc_tilt_coef``.
    """
    if INTERACTION_MODE == "contanc":
        m = CONTANC_M_S[study_1based - 1]
        return m, m
    if INTERACTION_MODE == "contanc_main":
        return +CONTANC_MAIN_DELTA, -CONTANC_MAIN_DELTA
    return 0.0, 0.0


# -----------------------------
# Covariate/Outcome simulation helpers
# -----------------------------

def _r_truncnorm(n: int, low: float, high: float, mean: float, sd: float, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """Sample n values from a truncated normal with bounds [low, high]."""
    rng = rng or np.random.default_rng()
    a, b = (low - mean) / sd, (high - mean) / sd
    return truncnorm.rvs(a, b, loc=mean, scale=sd, size=n, random_state=rng)


def simulate_population_df(
    X_std: np.ndarray,
    target_prev: float = 0.109,
    beta_age: float = 0.03,
    beta_gender: float = -0.2,
    # Genetic architecture
    h2: float = 0.0,
    num_causal: int = 0,
    scale_geffc_by_h2: bool = False,
    # Covariate distributions
    mode: str = 'homogeneous',  # 'homogeneous' or 'heterogeneous_A' or 'heterogeneous_B'
    # homogeneous age/gender params
    homog_age_low: float = 30.0,
    homog_age_high: float = 80.0,
    homog_gender_p: float = 0.5,
    # heterogeneous group A params (matches user's snippet)
    hetA_age_low: float = 30.0,
    hetA_age_high: float = 60.0,
    hetA_age_mean: float = 40.0,
    hetA_age_sd: float = 10.0,
    hetA_gender_p: float = 0.7,
    # heterogeneous group B params (matches user's snippet)
    hetB_age_low: float = 55.0,
    hetB_age_high: float = 80.0,
    hetB_age_mean: float = 70.0,
    hetB_age_sd: float = 10.0,
    hetB_gender_p: float = 0.3,
    # labels
    study_label: int = 0,
    population_label: Optional[str] = None,
    # RNG
    rng: Optional[np.random.Generator] = None,
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    """
    Create a DataFrame with covariates (age, gender), binary outcome Y, and labels
    for a *single* population given a standardized genotype matrix X_std (M x N).

    Parameters mirror the user's heterogeneous/homogeneous specifications. If
    `mode='heterogeneous_A'`, use (30–60) truncated-normal age + gender p=0.7.
    If `mode='heterogeneous_B'`, use (55–80) truncated-normal age + gender p=0.3.
    If `mode='homogeneous'`, use uniform age[30,80] + gender p=0.5.

    Genetic effect (geffc) defaults to zero unless you set `num_causal>0`.

    Returns
    -------
    df : pandas.DataFrame with columns ['age', 'gender', 'Y', 'study', 'population'] (population optional)
    extras : dict with keys {'geffc', 'beta_snp', 'causal_idx', 'p'} for reproducibility/inspection
    """
    rng = rng or np.random.default_rng()

    M, N = X_std.shape

    # --- Genetic effects
    beta_snp = np.zeros(M)
    causal_idx = np.array([], dtype=int)
    if num_causal and num_causal > 0:
        causal_idx = rng.choice(M, size=num_causal, replace=False)
        # Following the user's commented pattern:
        # draw effects ~ N(0, sqrt(h2/M)) and zero-out non-causal
        sd_beta = np.sqrt(h2 / max(M, 1)) if h2 > 0 else 0.0
        beta_draw = rng.normal(loc=0.0, scale=sd_beta, size=M)
        beta_snp[causal_idx] = beta_draw[causal_idx]
    geffc = X_std.T @ beta_snp  # shape (N,)

    # Optional re-scaling of geffc to match h2 relative to covariate variance
    # (matches idea in the user's snippet using v_g = h2 * v_cov / (1 - h2)).
    # We will compute v_cov after generating covariates and scale then.

    # --- Covariates
    if mode == 'homogeneous':
        age = rng.uniform(homog_age_low, homog_age_high, size=N)
        gender = rng.binomial(1, homog_gender_p, size=N)
    elif mode == 'heterogeneous_A':
        age = _r_truncnorm(N, hetA_age_low, hetA_age_high, hetA_age_mean, hetA_age_sd, rng)
        gender = rng.binomial(1, hetA_gender_p, size=N)
    elif mode == 'heterogeneous_B':
        age = _r_truncnorm(N, hetB_age_low, hetB_age_high, hetB_age_mean, hetB_age_sd, rng)
        gender = rng.binomial(1, hetB_gender_p, size=N)
    else:
        raise ValueError("mode must be one of {'homogeneous','heterogeneous_A','heterogeneous_B'}")

    # Variance components for optional scaling
    v_cov = (beta_age ** 2) * np.var(age) + (beta_gender ** 2) * np.var(gender)
    if scale_geffc_by_h2 and h2 > 0 and np.std(geffc) > 0:
        v_g = h2 * v_cov / (1 - h2)
        geffc = (geffc - geffc.mean()) / np.std(geffc) * np.sqrt(v_g)

    # --- Intercept to hit target prevalence
    logit_target = np.log(target_prev / (1 - target_prev))
    mean_cov_effect = beta_age * np.mean(age) + beta_gender * np.mean(gender)
    conf_eff = logit_target - mean_cov_effect - np.mean(geffc)

    # --- Outcome
    log_odds = geffc + conf_eff + beta_age * age + beta_gender * gender
    p = 1.0 / (1.0 + np.exp(-log_odds))
    Y = rng.binomial(1, p)

    # --- DataFrame
    df = pd.DataFrame({
        'age': age,
        'gender': gender,
        'Y': Y,
        'study': study_label,
    })
    if population_label is not None:
        df['population'] = population_label

    extras = {
        'geffc': geffc,
        'beta_snp': beta_snp,
        'causal_idx': causal_idx,
        'p': p,
    }
    return df, extras


def simulate_two_populations(
    X_std_A: np.ndarray,
    X_std_B: np.ndarray,
    heterogeneity: bool = False,
    target_prev: float = 0.109,
    target_prev_A: Optional[float] = None,
    target_prev_B: Optional[float] = None,
    beta_age: float = 0.03,
    beta_gender: float = -0.2,
    h2: float = 0.0,
    num_causal: int = 0,
    scale_geffc_by_h2: bool = False,
    rng: Optional[np.random.Generator] = None,
):
    """
    Convenience wrapper to create df (pop A) and df_2 (pop B) with either
    homogeneous or heterogeneous covariate distributions as requested.

    - If `heterogeneity=False`: both pops use the homogeneous scheme
      (age Uniform[30,80], gender p=0.5).
    - If `heterogeneity=True`: pop A uses truncated-normal 30–60, gender p=0.7;
      pop B uses truncated-normal 55–80, gender p=0.3.

    Returns
    -------
    df_A, df_B : DataFrames with columns ['age','gender','Y','study','population']
    info : dict with keys {'A': extrasA, 'B': extrasB}
    """
    rng = rng or np.random.default_rng()

    if heterogeneity:
        mode_A = 'heterogeneous_A'
        mode_B = 'heterogeneous_B'
    else:
        mode_A = mode_B = 'homogeneous'

    # Per-population prevalence: defaults pull from POP_A_PREV / POP_B_PREV,
    # but `target_prev` is honoured for back-compat if the per-pop kwargs
    # are not supplied (callers can override either independently).
    prev_A = target_prev_A if target_prev_A is not None else POP_A_PREV
    prev_B = target_prev_B if target_prev_B is not None else POP_B_PREV

    df_A, extra_A = simulate_population_df(
        X_std_A,
        target_prev=prev_A,
        beta_age=beta_age,
        beta_gender=beta_gender,
        h2=h2,
        num_causal=num_causal,
        scale_geffc_by_h2=scale_geffc_by_h2,
        mode=mode_A,
        study_label=0,
        population_label='PopA',
        rng=rng,
    )

    df_B, extra_B = simulate_population_df(
        X_std_B,
        target_prev=prev_B,
        beta_age=beta_age,
        beta_gender=beta_gender,
        h2=h2,
        num_causal=num_causal,
        scale_geffc_by_h2=scale_geffc_by_h2,
        mode=mode_B,
        study_label=1,
        population_label='PopB',
        rng=rng,
    )

    return df_A, df_B, {'A': extra_A, 'B': extra_B}

# -----------------------------
# Subpopulation Sampling Helpers
# -----------------------------
def _stratified_sample(
    indices: np.ndarray,
    y_values: np.ndarray,
    n_total: int,
    case_frac: float,
    rng: np.random.Generator,
    anc: Optional[np.ndarray] = None,
    m_case: float = 0.0,
    m_ctrl: float = 0.0,
) -> np.ndarray:
    """
    Draw ``n_total`` subjects from ``indices`` with a target case fraction.

    Separates indices into cases (Y=1) and controls (Y=0), draws the
    required number of each, and returns the selected indices.  If there
    aren't enough cases or controls, draws as many as available and fills
    the remainder from the other group.

    MECHANISM C (modes 7/8) only, via ``anc``.  When ``anc`` is given it is a
    per-panel array of continuous ancestry coordinates indexed by the SAME row
    positions as ``indices``, and ``m_case`` / ``m_ctrl`` are TARGET MEAN
    SHIFTS (in sd units of ``anc``) for the case draw and the control draw.
    Each pool is put in a WEIGHTED random order by the Gumbel-top-k identity --
    adding an i.i.d. Gumbel(0,1) to ``c*anc`` and sorting descending is exactly
    successive sampling with probability proportional to ``exp(c*anc)`` -- with
    the coefficient ``c`` SOLVED against that pool's own size and composition
    by ``contanc_tilt_coef`` so that the two pools reach the same shift even
    though the case pool is far more depleted than the control pool.  Passing
    the coefficient directly instead was measured to leak a within-study
    residual of -0.10 into an arm the design requires to be clean.

    When ``anc`` is None (every other mode, including OFF) the function is
    byte-for-byte the original: same two ``rng.shuffle`` calls in the same
    order, so the RNG stream and therefore the whole simulation is unchanged.
    """
    case_mask = y_values == 1
    case_idx = indices[case_mask]
    ctrl_idx = indices[~case_mask]
    if anc is None:
        rng.shuffle(case_idx)
        rng.shuffle(ctrl_idx)
    else:
        n_cases_target = int(round(n_total * case_frac))
        for pool, m, k in ((case_idx, m_case, n_cases_target),
                           (ctrl_idx, m_ctrl, n_total - n_cases_target)):
            if pool.size == 0:
                continue
            a_pool = anc[pool]
            c, _ach = contanc_tilt_coef(a_pool, min(k, pool.size), m)
            g = rng.gumbel(size=pool.size)
            order = np.argsort(-(c * a_pool + g), kind="stable")
            pool[:] = pool[order]

    n_cases_want = int(round(n_total * case_frac))
    n_ctrls_want = n_total - n_cases_want

    # Cap at available, redistribute if needed
    n_cases = min(n_cases_want, len(case_idx))
    n_ctrls = min(n_ctrls_want, len(ctrl_idx))
    shortfall = n_total - n_cases - n_ctrls
    if shortfall > 0:
        if n_cases < n_cases_want:
            n_ctrls = min(n_ctrls + shortfall, len(ctrl_idx))
        else:
            n_cases = min(n_cases + shortfall, len(case_idx))

    return np.concatenate([case_idx[:n_cases], ctrl_idx[:n_ctrls]])


def sample_subpopulations(
    df: pd.DataFrame,
    df_2: pd.DataFrame,
    rng: Optional[np.random.Generator] = None,
    anc_A: Optional[np.ndarray] = None,
    anc_B: Optional[np.ndarray] = None,
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Partition both populations into N_COHORTS=6 non-overlapping cohorts
    with **ancestry-skewed enrolment** AND **differential ascertainment**,
    mimicking the ADSP multi-site design.

    Per-study enrolment (see COHORT_POP_*_SIZES, COHORT_CASE_FRACTIONS):

      Study  PopA  PopB  case%   role
        1    475   25    80%     Pop A memory clinic   (Pop-A-skewed)
        2    475   25    20%     Pop A biobank
        3    250   250   50%     multi-ancestry biobank (mixed)
        4    250   250   80%     multi-ancestry clinic  (mixed)
        5     25  475    20%     Pop B biobank
        6     25  475    50%     Pop B clinic           (Pop-B-skewed)

    Total: 1500 Pop A + 1500 Pop B = 3000 subjects per replicate.

    Design rationale
    ----------------
    **Study ↔ ancestry correlation** (NEW vs prior design)
    Studies 1, 2 are 95% Pop A; studies 5, 6 are 95% Pop B; studies 3, 4
    are 50/50.  So ``study`` is partly predictive of ancestry — the same
    pattern as real ADSP, where seq-centre alone predicts PC1 well
    because some centres recruited disproportionately from one ancestry.
    This means methods that adjust ONLY for ``study`` (``combined_study``)
    will still leave residual stratification uncontrolled in some
    scenarios — but, unlike the previous balanced design, the residual
    is now real and detectable.

    **Differential case/control ascertainment (study-Y confound)**
    Six distinct case fractions (case-enriched / control-enriched /
    balanced × 2) preserve the ascertainment confound from the prior
    design.  All fractions are in (0, 1) so logistic regression remains
    identifiable.

    **Previous designs and why they were replaced**
    Design v1 (mono-ancestry studies):
      Studies 1/3/4 = Pop A only; 2/5/6 = Pop B only → ``study`` was a
      perfect proxy for ancestry.  Study 1 was case-only, study 2 was
      control-only → complete separation in logistic regression.
    Design v2 (multi-ancestry, random partition):
      Correctly made each study multi-ancestry, but assigned subjects
      to studies at random without regard to Y.  ``study`` was therefore
      independent of Y → ``combined_study`` faced no confounding and
      appeared trivially well-calibrated.

    Parameters
    ----------
    df   : DataFrame for Pop A (N_A rows).  Must have columns ``Y`` and
           ``study``, with a default integer RangeIndex.
    df_2 : DataFrame for Pop B (N_B rows).  Same requirements.
    rng  : numpy Generator; defaults to a fixed seed for reproducibility.

    Returns
    -------
    df_combined        : DataFrame of length <= 2 x N_COHORTS x COHORT_SIZE,
                         Pop A rows first then Pop B rows, with columns
                         ``study`` in {1, 2, 3} and ``population`` in
                         {PopA, PopB}.  Index is reset to 0-based.
    combined_indices   : 1-D int array of selected Pop A row positions.
    combined_indices_2 : Same for Pop B.
    """
    rng = rng or np.random.default_rng(RNG_SEED)

    assert len(COHORT_CASE_FRACTIONS) == N_COHORTS, (
        f"COHORT_CASE_FRACTIONS has {len(COHORT_CASE_FRACTIONS)} entries "
        f"but N_COHORTS={N_COHORTS}"
    )
    assert len(COHORT_POP_A_SIZES) == N_COHORTS and len(COHORT_POP_B_SIZES) == N_COHORTS, (
        f"COHORT_POP_*_SIZES must have N_COHORTS={N_COHORTS} entries"
    )

    # MECHANISM C switch.  False unless the mode is 7/8 AND the caller supplied
    # both continuous ancestry coordinates; when False nothing below changes.
    _contanc_on = (INTERACTION_MODE in _CONTANC_MODES
                   and anc_A is not None and anc_B is not None)
    if INTERACTION_MODE in _CONTANC_MODES and not _contanc_on:
        raise ValueError(
            f"mode {INTERACTION_MODE!r} requires anc_A and anc_B (the per-panel "
            f"continuous ancestry coordinates); got "
            f"{None if anc_A is None else anc_A.shape} / "
            f"{None if anc_B is None else anc_B.shape}.  Silently falling back "
            f"to the OFF sampler would mislabel an OFF run as a mode-7/8 run."
        )
    if _contanc_on:
        if len(anc_A) != len(df) or len(anc_B) != len(df_2):
            raise ValueError(
                f"ancestry coordinates must be aligned to the FULL panels: got "
                f"{len(anc_A)}/{len(anc_B)} for panels of {len(df)}/{len(df_2)}"
            )

    groups_A: List[np.ndarray] = []
    groups_B: List[np.ndarray] = []

    # Track which subjects have been assigned so far (no overlap across cohorts).
    used_A = np.zeros(len(df), dtype=bool)
    used_B = np.zeros(len(df_2), dtype=bool)

    for k in range(N_COHORTS):
        case_frac_A = COHORT_CASE_FRACTIONS_A[k]
        case_frac_B = COHORT_CASE_FRACTIONS_B[k]

        # Available subjects (not yet assigned to a previous cohort)
        avail_A = np.where(~used_A)[0]
        avail_B = np.where(~used_B)[0]

        n_per_A = min(COHORT_POP_A_SIZES[k], len(avail_A))
        n_per_B = min(COHORT_POP_B_SIZES[k], len(avail_B))

        # MECHANISM C (modes 7/8) only.  `_contanc_on` is False in every other
        # mode, so both calls below take their original 5-argument form and the
        # RNG stream is untouched.
        _m_case, _m_ctrl = (contanc_study_tilts(k + 1) if _contanc_on else (0.0, 0.0))
        _kw = dict(anc=anc_A, m_case=_m_case, m_ctrl=_m_ctrl) if _contanc_on else {}
        # Stratified draw: per-pop case fractions within the same study.
        sel_A = _stratified_sample(
            avail_A, df['Y'].values[avail_A], n_per_A, case_frac_A, rng, **_kw
        )
        if _contanc_on:
            _kw['anc'] = anc_B
        sel_B = _stratified_sample(
            avail_B, df_2['Y'].values[avail_B], n_per_B, case_frac_B, rng, **_kw
        )

        used_A[sel_A] = True
        used_B[sel_B] = True

        # Assign study label (1-indexed)
        df.loc[sel_A, 'study'] = k + 1
        df_2.loc[sel_B, 'study'] = k + 1

        groups_A.append(sel_A)
        groups_B.append(sel_B)

        # Log per-cohort composition
        cfA = df['Y'].values[sel_A].mean() if len(sel_A) > 0 else 0
        cfB = df_2['Y'].values[sel_B].mean() if len(sel_B) > 0 else 0
        print(
            f"  Study {k+1}: target_case_frac_A={case_frac_A:.0%} target_case_frac_B={case_frac_B:.0%}, "
            f"PopA={len(sel_A)} (case_frac={cfA:.1%}), "
            f"PopB={len(sel_B)} (case_frac={cfB:.1%})"
        )

    combined_indices = np.concatenate(groups_A)
    combined_indices_2 = np.concatenate(groups_B)

    n_total = len(combined_indices) + len(combined_indices_2)
    print(f"  Total subjects: {n_total}")

    df_combined = pd.concat(
        [df.iloc[combined_indices], df_2.iloc[combined_indices_2]],
        ignore_index=True,
    )
    return df_combined, combined_indices, combined_indices_2

# -----------------------------
# GMM Helper
# -----------------------------
def assign_population_labels(
    df: pd.DataFrame,
    n_components: int = 2,
    random_state: Optional[int] = None,
) -> pd.DataFrame:
    """
    Assign population labels using a GMM fitted on the first 20 PCs
    augmented with one-hot encoded ``center`` columns.  This lets the
    clustering account for sequencing-center batch structure in addition
    to genetic ancestry.

    Adds ``latent_class`` and ``assigned_pop`` columns to *df*.
    """
    pc_cols = [f'PC{i+1}' for i in range(N_PCS)]
    pcs = df[pc_cols].values

    center_dummies = pd.get_dummies(df['center'], prefix='ctr', dtype=float).values
    features = np.hstack([pcs, center_dummies])

    gmm = GaussianMixture(
        n_components=n_components,
        covariance_type='full',
        n_init=20,
        random_state=random_state,
    ).fit(features)

    latent_labels = gmm.predict(features)

    pop_labels = df['population'].values
    mapping = {}
    for k in range(n_components):
        majority_pop = pd.Series(pop_labels[latent_labels == k]).value_counts().idxmax()
        mapping[k] = majority_pop

    assigned = pd.Series(latent_labels).map(mapping).values

    acc = accuracy_score(pop_labels, assigned)
    cm = confusion_matrix(pop_labels, assigned, labels=np.unique(pop_labels))

    print(f"\nCenter-aware GMM assignment accuracy: {acc:.4f}")
    print("Confusion matrix (rows = true, cols = assigned):\n", cm)

    df['latent_class'] = latent_labels
    df['assigned_pop'] = assigned
    return df

# -----------------------------
# Batch-effect injection
# -----------------------------

def inject_batch_effects(
    X0: np.ndarray,
    center_labels: np.ndarray,
    frac_affected: float = 0.05,
    effect_size: float = 0.3,
    clip_dosage: bool = False,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """
    Simulate heterogeneous sequencing-center batch effects.

    Each center receives a *distinct* effect magnitude (``effect_size``
    scaled by a per-center multiplier spread across [0.5, 2.0]) and its
    own independently sampled set of affected SNPs whose fraction is
    also center-specific ([0.5x, 1.5x] of ``frac_affected``).  This
    guarantees statistically significant cross-center heterogeneity.

    Parameters
    ----------
    clip_dosage : bool
        If True, clip values to [0, 2] (appropriate for raw dosage data).
        Set False when operating on standardized genotypes.

    Returns a **copy** of X0 with artifacts applied (original is not mutated).
    """
    rng = rng or np.random.default_rng()
    M, N = X0.shape
    X0_obs = X0.astype(np.float64).copy()

    centers = np.unique(center_labels)
    n_centers = len(centers)

    es_multipliers = np.linspace(0.5, 2.0, n_centers)
    rng.shuffle(es_multipliers)
    frac_multipliers = np.linspace(0.5, 1.5, n_centers)
    rng.shuffle(frac_multipliers)

    for i, center in enumerate(centers):
        col_mask = center_labels == center
        es = effect_size * es_multipliers[i]
        n_affected = max(1, int(M * frac_affected * frac_multipliers[i]))
        snp_idx = rng.choice(M, size=n_affected, replace=False)
        directions = rng.choice([-1, 1], size=n_affected)
        for s, d in zip(snp_idx, directions):
            X0_obs[s, col_mask] += d * es

    if clip_dosage:
        X0_obs = np.clip(X0_obs, 0.0, 2.0)
    return X0_obs


def verify_batch_heterogeneity(
    X0_batch: np.ndarray,
    center_labels: np.ndarray,
    n_test_snps: int = 500,
    alpha: float = 0.05,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """
    Verify that injected batch effects produce statistically significant
    heterogeneity across sequencing centers via per-SNP one-way ANOVA
    (dosage ~ center).

    Asserts that > 10 % of sampled SNPs show a significant center effect,
    well above the 5 % nominal rate expected under the null.
    """
    rng = rng or np.random.default_rng()
    M, _N = X0_batch.shape
    centers = np.unique(center_labels)

    test_idx = rng.choice(M, size=min(n_test_snps, M), replace=False)
    sig_count = 0
    pvals: List[float] = []

    for snp in test_idx:
        groups = [X0_batch[snp, center_labels == c] for c in centers]
        groups = [g for g in groups if len(g) >= 2]
        if len(groups) < 2:
            continue
        _stat, p = f_oneway(*groups)
        pvals.append(float(p))
        if p < alpha:
            sig_count += 1

    n_tested = len(pvals)
    frac_sig = sig_count / max(n_tested, 1)
    median_p = float(np.median(pvals)) if pvals else 1.0

    print(f"  Batch heterogeneity: {sig_count}/{n_tested} SNPs "
          f"({frac_sig:.1%}) show significant center effects (ANOVA p<{alpha})")
    print(f"  Median ANOVA p-value: {median_p:.2e}")

    assert frac_sig > 0.10, (
        f"Batch effects insufficiently heterogeneous: only {frac_sig:.1%} "
        f"of tested SNPs significant — increase effect_size or frac_affected."
    )

    return {'frac_significant': frac_sig, 'median_pval': median_p, 'n_tested': n_tested}


# -----------------------------------------------------------------------------
# MECHANISM B — study-resolution batch component (INTERACTION_MODE 3 / 4 only)
# -----------------------------------------------------------------------------
# Everything below is inert unless INTERACTION_MODE is "studybatch" or
# "pseudostudy".  It never touches the caller's RNG in any other mode, so the
# off / center / study code paths remain bit-for-bit identical.
#
# Magnitude, in the units the code already uses: the study component is drawn
# with the SAME frac_affected (0.05) and effect_size (0.30) defaults as the
# centre component, and with the same multiplier ladders
# (linspace(0.5, 2.0, n_groups) on effect size, linspace(0.5, 1.5, n_groups) on
# affected fraction).  Because there are 6 studies and 3 centres, the ladders
# are finer but span the same range, and the measured added "mass"
# sum_g n_affected_g * es_g^2 * n_subjects_g is ~0.81x the centre component's.
# So mode 3 roughly DOUBLES the total batch artefact in the data — which is
# precisely why mode 4 exists.
STUDY_BATCH_FRAC = 0.05
STUDY_BATCH_SIZE = 0.30
# Dedicated RNG streams so modes 3 and 4 draw the IDENTICAL artefact (same SNP
# sets, same signs, same multipliers) and differ only in WHICH subjects receive
# it.  Offsets from the caller's seed, so the artefact still varies replicate
# to replicate.
STUDY_BATCH_ART_SEED = 1_000_003
STUDY_BATCH_PART_SEED = 2_000_003


def inject_group_batch_effects(
    X0: np.ndarray,
    group_labels: np.ndarray,
    frac_affected: float,
    effect_size: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, float]:
    """``inject_batch_effects`` generalised to an ARBITRARY grouping.

    Identical logic to ``inject_batch_effects`` (same multiplier ladders, same
    shuffle, same per-group independent SNP draw, same +/- directions) so the
    injected artefact is quantified in exactly the units the existing centre
    component uses.  Returns ``(X_copy, mass)`` where

        mass = sum_g  n_affected_g * es_g^2 * n_subjects_g

    is the total injected sum of squares, comparable across groupings.
    """
    M, _N = X0.shape
    X = X0.astype(np.float64).copy()
    groups = np.unique(group_labels)
    n_groups = len(groups)

    es_multipliers = np.linspace(0.5, 2.0, n_groups)
    rng.shuffle(es_multipliers)
    frac_multipliers = np.linspace(0.5, 1.5, n_groups)
    rng.shuffle(frac_multipliers)

    mass = 0.0
    for i, g in enumerate(groups):
        col_mask = group_labels == g
        es = effect_size * es_multipliers[i]
        n_affected = max(1, int(M * frac_affected * frac_multipliers[i]))
        snp_idx = rng.choice(M, size=n_affected, replace=False)
        directions = rng.choice([-1, 1], size=n_affected)
        for s, d in zip(snp_idx, directions):
            X[s, col_mask] += d * es
        mass += n_affected * es * es * int(col_mask.sum())
    return X, float(mass)


def pseudo_study_labels(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Permute ``study`` WITHIN each ``center``.

    Preserves the study -> centre map, the per-centre study-size composition,
    and therefore the group sizes the artefact is applied to.  Destroys the
    study <-> Y association (and the study <-> ancestry association, which is
    nil within a centre by design).
    """
    lab = df['study'].to_numpy().copy()
    ctr = df['center'].to_numpy()
    for c in np.unique(ctr):
        m = np.flatnonzero(ctr == c)
        lab[m] = lab[rng.permutation(m)]
    return lab


def maybe_inject_study_batch(
    X0: np.ndarray,
    df: pd.DataFrame,
    seed: int,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """Apply the Mechanism-B study batch component if the mode calls for it.

    NO-OP in modes off / center / study: returns ``X0`` unchanged (the same
    object, not a copy) and a diagnostics dict of zeros, and consumes nothing
    from any shared RNG.  Call it AFTER ``inject_batch_effects`` and BEFORE
    ``verify_batch_heterogeneity`` so the existing >10%-of-SNPs centre ANOVA
    assertion is evaluated on the final genotype matrix.
    """
    if INTERACTION_MODE not in _GENO_CHANNEL_MODES:
        return X0, {'study_batch_mass': 0.0, 'study_batch_applied': 0.0}

    rng_art = np.random.default_rng(int(seed) + STUDY_BATCH_ART_SEED)
    if INTERACTION_MODE == 'studybatch':
        labels = df['study'].to_numpy()
    else:  # 'pseudostudy'
        rng_part = np.random.default_rng(int(seed) + STUDY_BATCH_PART_SEED)
        labels = pseudo_study_labels(df, rng_part)

    X_out, mass = inject_group_batch_effects(
        X0, labels, STUDY_BATCH_FRAC, STUDY_BATCH_SIZE, rng_art
    )
    print(f"  [mechB] mode={INTERACTION_MODE} study-batch injected, "
          f"mass={mass:.4e}, n_groups={len(np.unique(labels))}")
    return X_out, {'study_batch_mass': mass, 'study_batch_applied': 1.0}


# Six studies → three sequencing centres, paired by ancestry profile so
# that each centre is dominated by one ancestry (mimicking ADSP, where
# Baylor/WashU/Broad/MEDGENOME differ markedly in ancestry composition):
#
#   Centre 1: studies 1, 2  (95% Pop A)        — Pop-A-dominant centre
#   Centre 2: studies 3, 4  (50% Pop A)        — multi-ancestry centre
#   Centre 3: studies 5, 6  (95% Pop B)        — Pop-B-dominant centre
#
# CENTER_NOISE_PROB controls how often a subject is reassigned to a
# non-primary centre.  At 0.002, study → centre concordance is ~99.8%,
# matching the ADSP empirical observation that study labels predict
# sequencing-centre essentially perfectly while leaving just enough
# noise to break complete collinearity in regression.
PRIMARY_STUDY_TO_CENTER = {1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3}
STUDY_TO_CENTER = PRIMARY_STUDY_TO_CENTER  # back-compat alias
# The interaction-design assertions near COHORT_CASE_FRACTIONS_A run before this
# definition exists, so they use a local copy of this map.  If the two ever
# diverge, the study-mode "every centre's mean gap is zero" guarantee silently
# stops holding, so pin them together here.
assert PRIMARY_STUDY_TO_CENTER == _STUDY_TO_CENTER_DESIGN, (
    f"PRIMARY_STUDY_TO_CENTER {PRIMARY_STUDY_TO_CENTER} no longer matches the "
    f"map the interaction design was verified against ({_STUDY_TO_CENTER_DESIGN}); "
    f"the study-mode centre-balance guarantee must be re-derived"
)
N_CENTERS = len(set(PRIMARY_STUDY_TO_CENTER.values()))
# Calibrated to give Study→Center ≈ 0.86 on out-of-sample multinomial
# logistic, matching the ADSP empirical Cohort-only prediction (85.9%).
# Each non-primary subject lands at one of the (N_CENTERS - 1) other
# centers uniformly at random.
CENTER_NOISE_PROB = 0.14


def _assign_centers(df: pd.DataFrame) -> pd.DataFrame:
    """Derive a ``center`` column from existing ``study`` labels.

    Maps each study to its primary centre via PRIMARY_STUDY_TO_CENTER
    (6 studies → 3 centres, two studies per centre).  A small fraction
    (CENTER_NOISE_PROB) of subjects are reassigned to a different centre
    so that ``study → center`` is ~99.8% but not perfectly collinear,
    matching the ADSP pattern (PC alone predicts centre at ~85%; study
    predicts centre at ~99.8%).

    Independent batch effects are injected per centre by
    ``inject_batch_effects``; each centre still processes subjects from
    both populations (mixed by design within Pop-A-dominant /
    Pop-B-dominant centres) so batch and ancestry are not perfectly
    collinear.
    """
    rng = np.random.default_rng(RNG_SEED + 1)
    primary = df['study'].map(PRIMARY_STUDY_TO_CENTER).fillna(0).astype(int).to_numpy()
    centers_arr = primary.copy()

    # Force exactly k subjects to be re-mapped so the noise is deterministic
    # (Bernoulli per subject is high-variance at small probabilities).
    n = len(centers_arr)
    n_flip = max(1, int(round(n * CENTER_NOISE_PROB)))
    flip_idx = rng.choice(n, size=n_flip, replace=False)
    all_centers = np.array(sorted(set(PRIMARY_STUDY_TO_CENTER.values())))
    for i in flip_idx:
        other = all_centers[all_centers != centers_arr[i]]
        centers_arr[i] = rng.choice(other)

    df['center'] = centers_arr
    return df

# -----------------------------
# Utilities
# -----------------------------

def _clean_pvals(p: pd.Series) -> np.ndarray:
    """Return valid p-values in (0,1]."""
    p = pd.to_numeric(p, errors='coerce').astype(float)
    p = p[(p > 0) & (p <= 1) & np.isfinite(p)]
    return p.values


def qq_plot_on_axes(ax: plt.Axes, pvals: np.ndarray, title: str):
    """Plot a single QQ curve onto ax (expected vs observed -log10 p)."""
    if pvals.size == 0:
        ax.text(0.5, 0.5, 'No valid p-values', ha='center', va='center', fontsize=10)
        ax.set_axis_off()
        return
    obs = -np.log10(np.sort(pvals))
    m = pvals.size
    exp = -np.log10((np.arange(1, m + 1) - 0.5) / m)
    maxv = max(exp.max(), obs.max())
    ax.plot(exp, obs, marker='.', linestyle='none', markersize=3)
    ax.plot([0, maxv], [0, maxv], linewidth=1)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('Expected -log10(p)')
    ax.set_ylabel('Observed -log10(p)')


def qq_grid(pvals_by_method: Dict[str, np.ndarray], fig_title: str, outfile: Optional[str] = None):
    """Make a 2x3 grid of QQ plots given a dict of method -> pvals."""
    method_order = [
        'combined_study',
        'combined_pcs',
        'combined_population',
        'combined_study_pcs',
        'combined_population_pcs',
        'matched_population'
    ]
    pretty_names = {
        'combined_study': '1) Combined + Study (+Age,Sex)',
        'combined_pcs': '2) Combined + PCs (+Age,Sex)',
        'combined_population': '3) Combined + Population (+Age,Sex)',
        'combined_study_pcs': '4) Combined + Study + PCs (+Age,Sex)',
        'combined_population_pcs': '5) Combined + Pop + PCs (+Age,Sex)',
        'matched_population': '6) Matched + Population (+Age,Sex)'
    }

    fig, axes = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()

    for i, key in enumerate(method_order):
        arr = pvals_by_method.get(key, np.array([]))
        qq_plot_on_axes(axes[i], arr, pretty_names[key])

    fig.suptitle(fig_title, fontsize=13)
    if outfile:
        fig.savefig(outfile, dpi=300)
    return fig

# -----------------------------
# Association engines (generic + wrappers)
# -----------------------------

def _prepare_covariate_matrix(df: pd.DataFrame, covariate_cols: List[str]) -> np.ndarray:
    """Build a numeric design matrix from covariates, with dummies and intercept."""
    Xc = df[covariate_cols].copy()
    for c in list(Xc.columns):
        if Xc[c].dtype.kind not in 'biufc':
            Xc[c] = Xc[c].astype('category')
        elif Xc[c].dtype.kind in 'iu' and Xc[c].nunique() > 2:
            Xc[c] = Xc[c].astype('category')
    Xc = pd.get_dummies(Xc, drop_first=True, dtype=float)
    Xc = sm.add_constant(Xc, has_constant='add')
    return Xc.values.astype(np.float64)


def assoc_Y_on_SNP_plus_covars(X0: np.ndarray, df: pd.DataFrame, covariate_cols: List[str]) -> pd.Series:
    """
    Vectorized logistic score (Rao) test: Y ~ covariates + SNP, H0: β_SNP = 0.

    Algorithm
    ---------
    1. Build the covariate design matrix W (categorical columns expanded to
       dummies, intercept appended).
    2. Fit the *null* logistic model  Y ~ W  once via statsmodels GLM/IRLS.
    3. Compute raw residuals  e = Y − μ  and working weights  w = μ(1−μ).
    4. Project out the covariate space from each SNP vector using a
       rank-revealing weighted QR decomposition of W (handles collinear
       dummies robustly).
    5. For each SNP g:
         score numerator   U  = g · e
         score variance    V  = ‖g_w‖² − ‖Q^T g_w‖²   (g_w = g·√w)
         test statistic    T  = U²/V  ~  χ²(1) under H0

    Convergence note
    ----------------
    If the null GLM does not converge (e.g. due to near-complete separation
    in a covariate), a warning is issued.  Complete separation—where a
    covariate perfectly predicts Y for a subset of subjects—drives the
    corresponding coefficient to ±∞ and the fitted probabilities for that
    subset to 0 or 1.  This makes w → 0 for those subjects, so they
    contribute near-zero weight to U and V and are effectively silently
    excluded from the score test.  To avoid this, ensure that no covariate
    (or covariate combination) perfectly separates cases from controls
    within the dataset passed to this function.

    Parameters
    ----------
    X0            : (M, N) genotype matrix (standardised dosages).
    df            : DataFrame with N rows aligned to X0 columns.
    covariate_cols: column names to include as covariates (categoricals
                    are dummy-encoded; an intercept is always added).

    Returns
    -------
    pd.Series of length M with p-values (NaN for degenerate SNPs).
    """
    import warnings

    y = df['Y'].values.astype(np.float64)
    W = _prepare_covariate_matrix(df, covariate_cols)

    M, n = X0.shape
    if n != len(df):
        raise ValueError(f"X0 columns (n={n}) must match df rows (N={len(df)}). Align first!")

    null_model = sm.GLM(y, W, family=sm.families.Binomial()).fit(disp=False)
    if not null_model.converged:
        warnings.warn(
            "Null GLM did not converge — score statistics may be unreliable. "
            "Check for complete or near-complete separation in the covariates: "
            f"{covariate_cols}",
            RuntimeWarning,
            stacklevel=2,
        )
    mu = null_model.predict()
    e = y - mu                              # (n,) raw residuals

    w = mu * (1.0 - mu)                     # (n,) working weights
    sqrt_w = np.sqrt(w)                     # (n,)

    W_tilde = W * sqrt_w[:, None]           # (n, p) weighted covariates
    Q, R, _ = scipy_qr(W_tilde, pivoting=True)
    rank = int(np.sum(np.abs(np.diag(R)) > 1e-8 * np.abs(R[0, 0])))
    Q = Q[:, :rank]

    G = X0.astype(np.float64)               # (M, n)

    U = G @ e                               # (M,) score numerator

    G_tilde = G * sqrt_w[None, :]           # (M, n) weighted SNPs
    GQ = G_tilde @ Q                        # (M, rank)
    V = np.sum(G_tilde * G_tilde, axis=1) - np.sum(GQ * GQ, axis=1)
    V[V < 1e-12] = np.nan

    score_stat = (U * U) / V               # chi-squared(1)
    pvals = chi2.sf(score_stat, df=1)
    pvals[~np.isfinite(pvals)] = np.nan

    return pd.Series(pvals, name='p_value')


# Method wrappers --------------------------------------------------------------

def pvals_combined_study(df: pd.DataFrame, X0: np.ndarray) -> pd.Series:
    return assoc_Y_on_SNP_plus_covars(X0, df, covariate_cols=['study', 'center', 'age', 'gender'])


def _pc_cols(df: pd.DataFrame) -> List[str]:
    pcs = [c for c in df.columns if c.upper().startswith('PC')]
    pcs_sorted = sorted(pcs, key=lambda x: int(''.join([ch for ch in x if ch.isdigit()]) or 0))
    return pcs_sorted


def pvals_combined_pcs(df: pd.DataFrame, X0: np.ndarray) -> pd.Series:
    pcs = _pc_cols(df)
    if not pcs:
        raise ValueError('No PC columns found (expected columns named PC1, PC2, ...).')
    covs = ['center', 'age', 'gender'] + pcs
    return assoc_Y_on_SNP_plus_covars(X0, df, covariate_cols=covs)


def pvals_combined_population(df: pd.DataFrame, X0: np.ndarray) -> pd.Series:
    return assoc_Y_on_SNP_plus_covars(X0, df, covariate_cols=['assigned_pop', 'center', 'age', 'gender'])


def pvals_combined_study_pcs(df: pd.DataFrame, X0: np.ndarray) -> pd.Series:
    pcs = _pc_cols(df)
    if not pcs:
        raise ValueError('No PC columns found (expected columns named PC1, PC2, ...).')
    covs = ['study', 'center', 'age', 'gender'] + pcs
    return assoc_Y_on_SNP_plus_covars(X0, df, covariate_cols=covs)


def pvals_combined_population_pcs(df: pd.DataFrame, X0: np.ndarray) -> pd.Series:
    pcs = _pc_cols(df)
    if not pcs:
        raise ValueError('No PC columns found (expected columns named PC1, PC2, ...).')
    covs = ['assigned_pop', 'center', 'age', 'gender'] + pcs
    return assoc_Y_on_SNP_plus_covars(X0, df, covariate_cols=covs)


# Matching-based method --------------------------------------------------------

def _fast_match(
    df: pd.DataFrame,
    group_cols: List[str],
    age_tolerance: int = 1,
) -> Tuple[List[int], List[int]]:
    """Greedy nearest-age matching with grouped index and binary search.

    For each case, looks up the control pool by *group_cols* in O(1) via a
    dict, then binary-searches the age-sorted pool for the nearest match
    within *age_tolerance*.  Matched controls are removed from the pool so
    each control is used at most once.  Overall O(N_cases * log(pool_size))
    instead of O(N_cases * N_controls).
    """
    cases = df[df['Y'] == 1]
    controls = df[df['Y'] == 0]

    pool: Dict[tuple, List[Tuple[float, int]]] = {}
    for idx, row in controls.iterrows():
        key = tuple(row[c] for c in group_cols)
        pool.setdefault(key, []).append((row['age'], idx))
    for v in pool.values():
        v.sort()

    matched_case: List[int] = []
    matched_ctrl: List[int] = []

    for case_idx, case in cases.iterrows():
        key = tuple(case[c] for c in group_cols)
        bucket = pool.get(key)
        if not bucket:
            continue
        target = case['age']
        pos = bisect_left(bucket, (target,))
        best_idx = None
        best_diff = float('inf')
        for candidate_pos in (pos - 1, pos):
            if 0 <= candidate_pos < len(bucket):
                diff = abs(bucket[candidate_pos][0] - target)
                if diff < best_diff:
                    best_diff = diff
                    best_idx = candidate_pos
        if best_idx is not None and best_diff <= age_tolerance:
            _, ctrl_idx = bucket.pop(best_idx)
            matched_case.append(case_idx)
            matched_ctrl.append(ctrl_idx)

    return matched_case, matched_ctrl

def pvals_matched_population(df: pd.DataFrame, X0: np.ndarray, age_tolerance: int = 1) -> pd.Series:
    cases, ctrls = match_cases_controls_with_center(df, age_tolerance=age_tolerance)
    if not cases or not ctrls:
        raise ValueError('Matching yielded empty sets. Check balance/criteria.')
    idx = cases + ctrls
    dfm = df.loc[idx].reset_index(drop=True)
    X0m = X0[:, idx]
    # Only adjust for age (inexactly matched, ±tolerance); exactly matched
    # variables (assigned_pop, center, gender) are already balanced and
    # including them wastes degrees of freedom / can cause near-separation.
    return assoc_Y_on_SNP_plus_covars(X0m, dfm, covariate_cols=['age'])


# Center-aware matching and GMM methods ----------------------------------------

def match_cases_controls_with_center(df: pd.DataFrame, age_tolerance: int = 1) -> Tuple[List[int], List[int]]:
    """Match on assigned_pop + center + gender + age using fast grouped lookup."""
    return _fast_match(df, group_cols=['gender', 'assigned_pop', 'center'], age_tolerance=age_tolerance)


# Study-resolution matching ----------------------------------------------------
# RIGA as published matches at CENTRE resolution.  Under
# SIM_INTERACTION_CONFOUNDING=2 the confounder lives at STUDY resolution and
# averages to exactly zero within every centre, so centre-level matching is
# blind to it.  These two arms exist to show that the failure is a property of
# the matching KEY, not of design-based correction as such.

def match_cases_controls_with_study(df: pd.DataFrame, age_tolerance: int = 1) -> Tuple[List[int], List[int]]:
    """Match on assigned_pop + study + gender + age using fast grouped lookup."""
    return _fast_match(df, group_cols=['gender', 'assigned_pop', 'study'],
                       age_tolerance=age_tolerance)


def match_cases_controls_with_study_center(df: pd.DataFrame, age_tolerance: int = 1) -> Tuple[List[int], List[int]]:
    """Match on assigned_pop + center + study + gender + age."""
    return _fast_match(df, group_cols=['gender', 'assigned_pop', 'center', 'study'],
                       age_tolerance=age_tolerance)


def pvals_matched_study(df: pd.DataFrame, X0: np.ndarray, age_tolerance: int = 1) -> pd.Series:
    """RIGA matched at STUDY resolution instead of centre resolution.

    Two caveats that must travel with any result from this arm:

    1.  `center` is NOT in the exact-match key, and centre is where
        `inject_batch_effects` puts the batch effects.  Study -> centre
        concordance is 1 - CENTER_NOISE_PROB, so this arm is exposed to batch
        confounding that `pvals_matched_population` is exactly balanced
        against.  `pvals_matched_study_center` is the companion that keeps both
        keys; use the pair to separate "matched at the right resolution" from
        "lost batch balance".
    2.  Study cells are half the size of centre cells, so this arm matches
        FEWER pairs (measured: ~0.61x of `matched_population` under the
        study-level manipulation).  Better calibration here is therefore
        partly an N effect and must be reported next to `n_matched_pairs`,
        never on its own.
    """
    cases, ctrls = match_cases_controls_with_study(df, age_tolerance=age_tolerance)
    if not cases or not ctrls:
        raise ValueError('Matching yielded empty sets. Check balance/criteria.')
    idx = cases + ctrls
    dfm = df.loc[idx].reset_index(drop=True)
    X0m = X0[:, idx]
    # Adjust only for age, exactly as pvals_matched_population does: the
    # exactly-matched variables are already balanced.
    return assoc_Y_on_SNP_plus_covars(X0m, dfm, covariate_cols=['age'])


def pvals_matched_study_center(df: pd.DataFrame, X0: np.ndarray, age_tolerance: int = 1) -> pd.Series:
    """Matched on assigned_pop + center + study (+ gender, age).

    The control for `pvals_matched_study`: it keeps the centre balance that
    protects against the injected batch effects while still resolving study, at
    the cost of the smallest matched sample of the three matched arms.
    """
    cases, ctrls = match_cases_controls_with_study_center(df, age_tolerance=age_tolerance)
    if not cases or not ctrls:
        raise ValueError('Matching yielded empty sets. Check balance/criteria.')
    idx = cases + ctrls
    dfm = df.loc[idx].reset_index(drop=True)
    X0m = X0[:, idx]
    return assoc_Y_on_SNP_plus_covars(X0m, dfm, covariate_cols=['age'])


# Matching PLUS PC covariates --------------------------------------------------
# The arms above adjust for age alone, on the argument that the exactly-matched
# variables are already balanced.  That argument holds only for variation at the
# RESOLUTION of the matching key.  A genotype gradient that varies WITHIN an
# assigned_pop label is not balanced by matching on the label, and the two arms
# below are the ones that say whether standard PC covariates mop it up.
#
# These are also the arms that correspond to the method this project actually
# deploys on ADSP: `lpa_pc` is LPA-cluster matching PLUS PC covariates, not
# matching alone.  Reporting a matching-only arm as "RIGA as published" against
# regression arms that all carry PC1-PC20 compares matching-without-PCs to
# regression-with-PCs, which is not the contrast the paper makes.
#
# The PCs are the FULL-SAMPLE PCs subset to the matched rows.  PCA is not refit
# on the matched subsample, because the deployed pipeline computes PCs once on
# the analysis sample and carries them into the matched model.

def _pvals_matched_with_pcs(df: pd.DataFrame, X0: np.ndarray,
                            matcher, age_tolerance: int = 1) -> pd.Series:
    cases, ctrls = matcher(df, age_tolerance=age_tolerance)
    if not cases or not ctrls:
        raise ValueError('Matching yielded empty sets. Check balance/criteria.')
    idx = cases + ctrls
    dfm = df.loc[idx].reset_index(drop=True)
    X0m = X0[:, idx]
    pcs = _pc_cols(df)
    if not pcs:
        raise ValueError('No PC columns found (expected columns named PC1, PC2, ...).')
    return assoc_Y_on_SNP_plus_covars(X0m, dfm, covariate_cols=['age'] + pcs)


def pvals_matched_population_pcs(df: pd.DataFrame, X0: np.ndarray,
                                 age_tolerance: int = 1) -> pd.Series:
    """RIGA's centre-resolution matching plus PC covariates (the deployed form)."""
    return _pvals_matched_with_pcs(df, X0, match_cases_controls_with_center,
                                   age_tolerance=age_tolerance)


def pvals_matched_study_pcs(df: pd.DataFrame, X0: np.ndarray,
                            age_tolerance: int = 1) -> pd.Series:
    """Study-resolution matching plus PC covariates."""
    return _pvals_matched_with_pcs(df, X0, match_cases_controls_with_study,
                                   age_tolerance=age_tolerance)

# -----------------------------
# Run all methods for one setting
# -----------------------------

def run_all_methods_for_setting(df: pd.DataFrame, X0: np.ndarray,
                                age_tolerance: int = 1) -> Dict[str, np.ndarray]:
    """Return dict of method_key -> pval ndarray for a single setting."""
    out: Dict[str, np.ndarray] = {}

    # Method 1
    p1 = pvals_combined_study(df, X0)
    out['combined_study'] = _clean_pvals(p1)

    # Method 2
    p2 = pvals_combined_pcs(df, X0)
    out['combined_pcs'] = _clean_pvals(p2)

    # Method 3
    p3 = pvals_combined_population(df, X0)
    out['combined_population'] = _clean_pvals(p3)

    # Method 4
    p4 = pvals_combined_study_pcs(df, X0)
    out['combined_study_pcs'] = _clean_pvals(p4)

    # Method 5
    p5 = pvals_combined_population_pcs(df, X0)
    out['combined_population_pcs'] = _clean_pvals(p5)

    # Method 6 (matched)
    p6 = pvals_matched_population(df, X0, age_tolerance=age_tolerance)
    out['matched_population'] = _clean_pvals(p6)

    # Methods 7-8 (matched at study resolution).  Guarded: their exact-match
    # cells are 2x and 4x finer than method 6's, so they can legitimately
    # produce an empty match on a small or badly balanced replicate and must
    # not be allowed to kill the whole setting.
    for _name, _fn in (('matched_study', pvals_matched_study),
                       ('matched_study_center', pvals_matched_study_center),
                       ('matched_population_pcs', pvals_matched_population_pcs),
                       ('matched_study_pcs', pvals_matched_study_pcs)):
        try:
            out[_name] = _clean_pvals(_fn(df, X0, age_tolerance=age_tolerance))
        except Exception as exc:                                   # noqa: BLE001
            print(f"  [{_name}] failed: {type(exc).__name__}: {exc}")
            out[_name] = np.array([])

    return out


# -----------------------------
# Driver: build figures for four settings
# -----------------------------

def make_all_qq_figures(settings: Dict[str, Tuple[pd.DataFrame, np.ndarray]],
                        age_tolerance: int = 1,
                        out_prefix: str = 'qq_') -> Dict[str, Dict[str, np.ndarray]]:
    """
    settings: dict mapping setting_key -> (df, X0)
      where setting_key is used in the output filename & title.
    Returns a nested dict of p-values by setting and method for any downstream use.
    """
    collected = {}
    for setting_key, (df, X0) in settings.items():
        print(f"Running six-method GWAS for setting: {setting_key}")
        pvals_by_method = run_all_methods_for_setting(df, X0, age_tolerance)
        collected[setting_key] = pvals_by_method
        title = (
            f"QQ plots — {setting_key.replace('_', ' ')}\n"
            "Methods: (1) Combined+Study, (2) Combined+PCs, (3) Combined+Pop,\n"
            "(4) Combined+Study+PCs, (5) Combined+Pop+PCs, (6) Matched+Pop"
        )
        outfile = f"{out_prefix}{setting_key}.png"
        qq_grid(pvals_by_method, fig_title=title, outfile=outfile)
        print(f"  → saved {outfile}")
    return collected


def _build_setting_dataset(
    X_std_A: np.ndarray,
    X_std_B: np.ndarray,
    heterogeneity: bool,
    n_components: int,
    seed: int,
    batch_effect: bool = True,
    batch_frac: float = 0.05,
    batch_size: float = 0.3,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Build one analysis-ready setting:
      1) simulate two populations (phenotypes from Xstd which carries LD structure)
      2) sample study subpopulations
      3) combine LD-structured genotypes (Xstd)
      4) inject heterogeneous sequencing-center batch effects
      5) verify batch heterogeneity (ANOVA across centers)
      6) compute PCs on concatenated Xstd (sklearn PCA centers internally)
      7) assign population labels (plain GMM + center-aware GMM)

    Both PCA and GWAS regression use Xstd (LD-structured genotypes).
    """
    rng = np.random.default_rng(seed)
    df_a, df_b, _ = simulate_two_populations(
        X_std_A=X_std_A,
        X_std_B=X_std_B,
        heterogeneity=heterogeneity,
        h2=0.0,
        num_causal=0,
        scale_geffc_by_h2=False,
        rng=rng,
    )
    # MECHANISM C (modes 7/8): the continuous ancestry coordinate must exist
    # BEFORE enrolment, because it is what enrolment is tilted on.  It is a
    # deterministic function of the panels (cached, fixed seed, no draw from
    # `rng`), so modes off/1-6 are bit-identical -- the two names stay None and
    # `sample_subpopulations` takes its original call form.
    _anc_A = _anc_B = None
    if INTERACTION_MODE in _CONTANC_MODES:
        _anc_A, _anc_B = contanc_scores(X_std_A, X_std_B)
    df_s, idx_a, idx_b = sample_subpopulations(
        df_a, df_b, rng=rng, anc_A=_anc_A, anc_B=_anc_B
    )
    if _anc_A is not None:
        # Carried through for the probe only.  It is never a covariate and is
        # not picked up by `_pc_cols` (which matches names starting "PC").
        df_s['anc_coord'] = np.concatenate([_anc_A[idx_a], _anc_B[idx_b]])
    Xstd_s = np.concatenate((X_std_A[:, idx_a], X_std_B[:, idx_b]), axis=1)

    if batch_effect:
        df_s = _assign_centers(df_s)
        center_labels = df_s['center'].values
        Xstd_s = inject_batch_effects(
            Xstd_s, center_labels,
            frac_affected=batch_frac,
            effect_size=batch_size,
            clip_dosage=False,
            rng=rng,
        )
        # Mechanism B (modes 3/4 only): no-op otherwise, and it does not touch
        # `rng`, so the centre-only path stays bit-identical.  Placed BEFORE the
        # ANOVA check so the existing centre assertion runs on the final matrix.
        Xstd_s, _mechb = maybe_inject_study_batch(Xstd_s, df_s, seed)
        verify_batch_heterogeneity(Xstd_s, center_labels, rng=rng)

    pca = PCA(n_components=N_PCS, svd_solver='randomized', random_state=seed)
    pcs = pca.fit_transform(Xstd_s.T)
    pcs_df = pd.DataFrame(pcs, columns=[f'PC{i+1}' for i in range(pcs.shape[1])])
    df_s = pd.merge(df_s, pcs_df, left_index=True, right_index=True)
    df_s = assign_population_labels(df_s, n_components=n_components, random_state=seed)

    return df_s, Xstd_s


def _scenario_configs(arrays: Dict[str, np.ndarray]) -> List[Dict[str, object]]:
    """
    Define simulation scenarios with labels/order used in the paper figure.
    All scenarios include heterogeneous batch effects whose cross-center
    heterogeneity is verified via one-way ANOVA (see verify_batch_heterogeneity).
    """
    return [
        {
            'key': '2pure_homo',
            'label': '2 pure populations, homogeneous covariates',
            'X_std_A': arrays['Xstd_afr'],
            'X_std_B': arrays['Xstd_eur'],
            'heterogeneity': False, 'n_components': 2,
            'batch_effect': True, 'batch_frac': 0.05, 'batch_size': 0.3,
        },
        {
            'key': '2pure_hetero',
            'label': '2 pure populations, heterogeneous covariates',
            'X_std_A': arrays['Xstd_afr'],
            'X_std_B': arrays['Xstd_eur'],
            'heterogeneity': True, 'n_components': 2,
            'batch_effect': True, 'batch_frac': 0.05, 'batch_size': 0.3,
        },
        {
            'key': '1pure_1admixed_homo',
            'label': '1 pure + 1 admixed, homogeneous covariates',
            'X_std_A': arrays['Xstd_afr'],
            'X_std_B': arrays['Xstd_ceu5yri5'],
            'heterogeneity': False, 'n_components': 2,
            'batch_effect': True, 'batch_frac': 0.05, 'batch_size': 0.3,
        },
        {
            'key': '1pure_1admixed_hetero',
            'label': '1 pure + 1 admixed, heterogeneous covariates',
            'X_std_A': arrays['Xstd_afr'],
            'X_std_B': arrays['Xstd_ceu5yri5'],
            'heterogeneity': True, 'n_components': 2,
            'batch_effect': True, 'batch_frac': 0.05, 'batch_size': 0.3,
        },
        # ── 2-admixed note ────────────────────────────────────────────────────
        # CEU5-YRI5 (50 % CEU + 50 % YRI) and CEU5-CHS5 (50 % CEU + 50 % CHS)
        # share half their ancestry (CEU) and are simulated from chr22 only.
        # Chr22 carries insufficient ancestry-informative markers to separate
        # the YRI vs CHS components in PCA.  As a result:
        #   • The two clouds overlap extensively in PC1/PC2 space.
        #   • GMM population assignment is near-chance (accuracy ≈ 50 %).
        #   • All six GWAS methods perform similarly — there is no detectable
        #     population structure to correct for with the available data.
        # This is an expected scientific finding, not a software defect.
        # ─────────────────────────────────────────────────────────────────────
        {
            'key': '2admixed_homo',
            'label': '2 admixed populations, homogeneous covariates',
            'X_std_A': arrays['Xstd_ceu5yri5'],
            'X_std_B': arrays['Xstd_ceu5chs5'],
            'heterogeneity': False, 'n_components': 2,
            'batch_effect': True, 'batch_frac': 0.05, 'batch_size': 0.3,
        },
        {
            'key': '2admixed_hetero',
            'label': '2 admixed populations, heterogeneous covariates',
            'X_std_A': arrays['Xstd_ceu5yri5'],
            'X_std_B': arrays['Xstd_ceu5chs5'],
            'heterogeneity': True, 'n_components': 2,
            'batch_effect': True, 'batch_frac': 0.05, 'batch_size': 0.3,
        },
    ]


def run_single_simulation_job(
    arrays: Dict[str, np.ndarray],
    replicate_id: int,
    seed: int,
    age_tolerance: int,
    out_dir: str,
):
    """
    Run one replicate (intended for one SLURM task) and save:
      - run_{replicate_id:03d}_metrics.csv
      - run_{replicate_id:03d}_pvals.npz
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    scenarios = _scenario_configs(arrays)
    metrics_rows: List[Dict[str, object]] = []
    pvals_archive: Dict[str, np.ndarray] = {}

    for s_idx, cfg in enumerate(scenarios):
        scenario_seed = seed + s_idx
        df_s, Xstd_s = _build_setting_dataset(
            X_std_A=cfg['X_std_A'],
            X_std_B=cfg['X_std_B'],
            heterogeneity=cfg['heterogeneity'],
            n_components=cfg['n_components'],
            seed=scenario_seed,
            batch_effect=cfg.get('batch_effect', False),
            batch_frac=cfg.get('batch_frac', 0.05),
            batch_size=cfg.get('batch_size', 0.3),
        )
        pvals_by_method = run_all_methods_for_setting(df_s, Xstd_s, age_tolerance=age_tolerance)
        for method, pvals in pvals_by_method.items():
            pvals_archive[f"{cfg['key']}__{method}"] = pvals
            metrics_rows.append({
                'replicate': replicate_id,
                'seed': seed,
                'scenario_key': cfg['key'],
                'scenario_label': cfg['label'],
                'method': method,
                'n_valid_pvals': int(len(pvals)),
                'lambda_gc': compute_lambda_gc(pvals),
                'type1_error_5e8': type1_error_rate(pvals, 5e-8),
                'type1_error_5e7': type1_error_rate(pvals, 5e-7),
                'type1_error_1e5': type1_error_rate(pvals, 1e-5),
            })

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_path = out_path / f"run_{replicate_id:03d}_metrics.csv"
    pvals_path = out_path / f"run_{replicate_id:03d}_pvals.npz"
    metrics_df.to_csv(metrics_path, index=False)
    np.savez_compressed(pvals_path, **pvals_archive)
    print(f"Saved run outputs: {metrics_path}, {pvals_path}")


def _load_default_arrays() -> Dict[str, np.ndarray]:
    """Load LD-structured standardized genotype arrays used by all SLURM tasks."""
    return {
        'Xstd_afr': np.load('${ADSP_ROOT}/ldsc_simu_std_afr.npy'),
        'Xstd_eur': np.load('${ADSP_ROOT}/ldsc_simu_std_eur.npy'),
        'Xstd_ceu5yri5': np.load('${ADSP_ROOT}/data/example_data/CEU5-YRI5_multichr_simu_std.npy'),
        'Xstd_ceu5chs5': np.load('${ADSP_ROOT}/data/example_data/CEU5-CHS5_multichr_simu_std.npy'),
        'Xstd_eur2': np.load('${ADSP_ROOT}/ldsc_simu_std.npy'),
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run one simulation replicate and write per-run output files."
    )
    parser.add_argument("--replicate-id", type=int, default=1)
    parser.add_argument("--seed", type=int, default=RNG_SEED)
    parser.add_argument("--age-tolerance", type=int, default=1)
    parser.add_argument("--out-dir", type=str, default="outputs/runs")
    return parser.parse_args()


# -----------------------------
# SETTINGS (FILL THESE IN)
# -----------------------------
if __name__ == '__main__':
    args = _parse_args()
    arrays = _load_default_arrays()
    run_single_simulation_job(
        arrays=arrays,
        replicate_id=args.replicate_id,
        seed=args.seed,
        age_tolerance=args.age_tolerance,
        out_dir=args.out_dir,
    )

