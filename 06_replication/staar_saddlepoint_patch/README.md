# ukb_rootcheck_20261019 — root-checked STAAR/STAARpipeline library and UK Biobank coding re-scan

Prepared 19 October 2026 (NOT submitted); **rev 2** after two reviews (section 10 lists the fixes). Everything in this
package is new; no existing file under `revision/scripts`, `revision/ukb`, `revision/tables`,
`revision/five_analysis_execution_20261016` or the default R library was modified. Submit with
`bash SUBMIT.sh [external job ids...]` (see "Submission").

## 1. Why

The installed STAAR 0.9.8 (xihaoli/STAAR commit `9db9dd50…`) and STAARpipeline 0.9.8 (`7e868f23…`) share a
saddlepoint solver whose Newton routine (`NR_Binary_SPA`) returns a single number with no convergence flag. The callers
(`STAAR_B_Binary_SPA`, `STAAR_B_SPA`, `STAAR_B_SPA_SMMAT`, `STAAR_B_SPA_SMMAT_sparse`; STAARpipeline's private copy in
`Individual_Score_Test_SPA.cpp`) fall back to bisection only when the resulting p is NA or exactly 1, so an unconverged
iterate that yields a finite tiny p passes through. Documented case (AUDIT_five_analyses_20261019.md §5): UK Biobank SNX22
synonymous, STAAR-B 2.8e-16, true value 0.052 (root 20.6 with score-equation residual 1.53 instead of root 4.21).
The bisection fallback had its own silent failure: when the endpoints ±100 did not bracket the root it returned x = 0,
which the callers reported as p = 1 (the component then silently dropped out of the Cauchy combination).

## 2. The patch (patches/*.patch, applied by build_patched_staar.sh)

Unified diffs against the untrusted upstream trees in `src_untrusted/` (GitHub archive tarballs by commit SHA; sha256 in
`build/BUILD_RECORD.txt`). Exported C++ signatures are unchanged, so `RcppExports.*` are untouched. The patches
re-create `build_work/` exactly (verified with `patch -p1` on fresh copies + `diff -r`).

STAAR (`patches/STAAR_saddlepoint_root_check.patch`, 10 files):
* `src/NR_Binary_SPA.cpp` — the original iteration is kept verbatim; afterwards the score-equation residual
  |K1(x) − q| is evaluated at the returned iterate. A root with |residual| ≤ 1e-6·max(1,|q|) is returned exactly as
  before (bit-identical p). A non-finite iterate, |x| > 1e10, or max_iter reached → `NA_REAL`. If the loop stopped on
  its own (absolute, 1.2e-4) criterion but the residual exceeds the relative tolerance, up to 20 further Newton steps
  polish the root from the same iterate; if that fails → `NA_REAL`.
* `src/Bisection_Binary_SPA.cpp` — non-finite endpoints → NA; endpoints of the same sign → the bracket is expanded
  outward by doubling (K1 − q is increasing in x), NA if no sign change after 60 doublings; bisection then runs until
  |K1(x) − q| ≤ 1e-6·max(1,|q|) (or the interval is exhausted at machine precision); the result is returned only if it
  passes that residual check, otherwise NA.
* `src/Saddle_Binary_SPA.cpp`, `src/Saddle_Binary_SPA_Bisection.cpp` — return NA when the root is NA (no formula
  evaluated on a NaN root).
* the four `STAAR_B_*SPA*.cpp` callers — the existing fallback (Newton NA or p = 1 → bisection) is kept; a component for
  which BOTH solvers failed is now kept as NA instead of 1 (`STAAR_Binary_SPA.R` already drops NA components from the
  Cauchy combination and leaves NA in the component column, so the failure is visible in the output); `res(i) > 1` is
  guarded against NaN. Each caller records its context for the diagnostic log and writes one `call_done` line (rev 2).
* new `src/rootcheck.cpp` — the tolerance, a K1 helper, and the diagnostic log. If the environment variable
  `STAAR_ROOTCHECK_LOG` names a file, one tab-separated line is appended per solver EVENT (Newton failure, Newton
  polish, bisection used/failed) and one `call_done` line per burden call; nothing is written for an ordinary verified
  Newton root. **17 columns (rev 2):** pkg, call_seq, n_variants, cMAC, n_weights, component, p_normal, tail, event, q,
  x_newton, resid_newton, x_final, resid_final, iter, **context, n_spa**. `component` is the 0-based weight column of the
  call (stored order: Burden(1,25), Burden(1,25)-CADD, …, Burden(1,1), Burden(1,1)-CADD, …; 26 columns for every
  coding category, disruptive_missense included); `p_normal` is the normal-approximation p that sent the component into
  SPA (NA for the unfiltered caller); `context` is the value of `STAAR_ROOTCHECK_CONTEXT` at the start of the call —
  the patched STAARpipeline sets it to `<gene>|<category>` — so every event is attributable to exactly one result row;
  `n_spa` (call_done lines) is the number of components of the call that entered the saddlepoint branch, the
  denominator of the root-failure census. Cost: a few file appends per burden call. Not thread-safe (the scan is 1 CPU).
* `DESCRIPTION` — a `RootcheckPatch:` field so that `packageDescription("STAAR")$RootcheckPatch` identifies the build.

STAARpipeline (`patches/STAARpipeline_saddlepoint_root_check.patch`, 8 files): the same solver changes inside the private copy
`src/Individual_Score_Test_SPA.cpp` (helpers file-local; log lines carry pkg = STAARpipeline, cMAC = the variant's
allele count, the same 17 columns); the `DESCRIPTION` field; and (rev 2) six R files — `coding.R` (the function behind
`Gene_Centric_Coding(category = "all_categories")`, which the scan uses) and the single-category functions `plof.R`,
`plof_ds.R`, `missense.R` (both of its calls), `disruptive_missense.R`, `synonymous.R` — gain
`Sys.setenv(STAAR_ROOTCHECK_CONTEXT = paste0(gene_name, "|<category>"))` immediately before each `STAAR_Binary_SPA`
call and an `on.exit(Sys.unsetenv(...))`. These lines touch no computation; the validation reruns prove it (every
event-free column bit-identical to the original scan).

Tolerance rationale: the relative error of a Lugannani–Rice p is ≈ |K1(x) − q| / |q| (w is stationary at the root), so
1e-6·max(1,|q|) verifies every p to about 1e-6 relative. Original Newton stops accepted |K1 − q| ≤ 1.2e-4 absolute;
roots that were accurate to between 1e-6 and 1.2e-4 are "polished" (logged as `newton_polished`) and their p values
move by at most ~1e-4 relative — these show up as `tiny_change` in the comparisons and are not failures.

Facts about the pipeline that the comparison relies on (read from `build_work/STAARpipeline/R/coding.R`):
* `coding()` calls `STAAR_Binary_SPA` once per category, **in the order plof_ds, plof, synonymous, missense,
  disruptive_missense** (the output list order is plof, plof_ds, missense, disruptive_missense, synonymous), always with
  the full annotation matrix (26 weight columns; disruptive_missense is NOT a 2-weight call — rev 1 had this wrong).
* The single disruptive_missense call feeds two rows: its own row and the missense row's `Burden(1,25)-Disruptive` /
  `Burden(1,1)-Disruptive` columns (weight columns 0 and 13 of that call, copied by coding.R; NA → 1 there), and the
  missense row's STAAR-B(1,25)/STAAR-B(1,1)/STAAR-B aggregates include those two copies.
* (#SNV, cMAC) is not unique within a chunk (plof == plof_ds whenever a gene has no disruptive missense variant; chunk
  676 has such pairs), which is why attribution uses the context and only VERIFIES (#SNV, cMAC).

## 3. Library

`build_patched_staar.sh` (run by `jobs/01_build.sbatch`; it was also run in the login shell on 19 Oct — rev 1 at 12:15, rev 2 at
13:23, see `build/BUILD_RECORD.txt` and `build_login_shell_rev2.log`; STAAR ~2.5 min + STAARpipeline ~2 min at -j2,
< 1 GB) installs both packages into the NEW library `${R_LIBS_ROOT}/staar_rootcheck_4.4`. Only the job
scripts of this package prepend it: `env_rootcheck.sh` sets
`R_LIBS_USER=${R_LIBS_ROOT}/staar_rootcheck_4.4:<default library>`; every other package resolves from the
default library as before. Installed versions: STAAR 0.9.8 + RootcheckPatch 2026-10-19 rev 2 (upstream
9db9dd504905b9f469146f670e5f6dbe3e08d01a), STAARpipeline 0.9.8 + RootcheckPatch 2026-10-19 rev 2 (upstream
7e868f23c85ba0713dd00c8512ad23cca7fbc8f9). The default library's RemoteSha fields are checked by build_patched_staar.sh.

## 4. Validation gate (jobs/02_validate.sbatch → VALIDATION.txt, gate/ROOTCHECK_GATE.json)

* (d) `scripts/unit_test.R <mode> <outdir> [null.RData]` — synthetic cases, run with the DEFAULT library
  (`original`: documents the defect) and the root-checked one (`patched`: the test proper):
  - **A** (solver level: one weight-16 variant, 900 weight-1 variants, heterogeneous μ, q = 11.56): the ORIGINAL Newton
    returns root 0.554 with residual 14.4 and p = 4.8e-219; the patched Newton returns NA, the patched bisection root
    has residual 9.8e-7 and p = 0.03166, equal to an independent `uniroot` + Lugannani–Rice reference to 1e-6.
  - **B** (caller level: fitted logistic null with a covariate, 5,000 subjects, 901 single-carrier variants, variant 1
    with weight 22.9 has one case and one control carrier): the ORIGINAL `STAAR_B_Binary_SPA` returns p = 1.76e-5 where
    the correct two-sided value is 0.0283 (the SNX22 pathology); the patched caller returns the bisection value, equal
    to the reference to 1e-6. **B4 (rev 2):** the same construction through the three other callers — `STAAR_B_SPA`
    (GLM, filtered), `STAAR_B_SPA_SMMAT` (dense P) and **`STAAR_B_SPA_SMMAT_sparse`, the caller the scan uses**
    (relatedness = TRUE, sparse_kins = TRUE, SPA_p_filter = TRUE) — with the null expressed as glmmkin quantities
    (Σ⁻¹ = diag(W)); all four must return the same bisection value to 1e-8, and the log must show one `call_done` line
    per caller with `n_spa` = 1 and `p_normal` < 0.05 for the filtered callers.
  - **C**: a target outside the support gives NA at the solver level. **G (rev 2):** a no-root construction that
    reaches every caller (an artificial residual vector of 3 for every subject makes the score exceed sup K1): all four
    callers must return **NA in the component column** (not 1, not a number), and the log must show the Newton failure
    and the bisection failure for each.
  - **D**: a well-conditioned case returns the Newton root with no event logged. **E**: STAARpipeline's single-variant
    SPA copy on the case-B construction.
  - **H (rev 2, when the null file is given — the validation job and the login run pass the reconstructed null):** the
    scan's own path, `STAAR_Binary_SPA()` with the REAL 182,238-subject sparse-kins null, SPA_p_filter = TRUE: 301
    variants, variant 1 (annotation PHRED 40 → weight ≈ 25) with one case and one control carrier, 300 singletons
    (PHRED 0.2 → weight ≈ 1.1), three of them in cases. Found by scanning constructions on the null on 19 Oct: both
    Newton tails of the annotation-weighted columns hit max_iter (reference two-sided p 0.0046). The patched caller must
    log the failure, fall back to bisection and equal the reference to 1e-6 for every SPA component; the `call_done`
    line must carry n_variants 301, cMAC 302, n_weights 4, n_spa = the number of SPA components. Forced NA: a copy of
    the null whose scaled residuals are 3 everywhere → every component NA in the component column (STAAR-B falls to 1
    by the unchanged R logic "all components NA → 1"; the audit counts NA components).
  **Login-shell run on 19 Oct with the rev-2 library: see section 11 and `results/validation_login/unit_test_rev2.log`.**
* (a) `scripts/validate.R` recomputes SNX22 synonymous, TREM2 missense and NECTIN2 missense from the exact
  masks cached by the 18 Oct targeted run (same null, G, annotations, disruptive subset, Cauchy combination) and,
  independently, through `Gene_Centric_Coding(category = …)` on the aGDS; both must agree with the references
  0.0518593609546751 / 2.10583109606644e-05 / 0.000714815559568427 to 1e-6 relative (components to 1e-5 relative against
  `components_guarded.csv`). The exact-mask half was run in the login shell (rev 1 and again with rev 2,
  `results/validation_login/refsets_exact_mask_login.{R,log}`, `refsets_rootlog.tsv`): SNX22 synonymous STAAR-B
  0.0518593585025056 (rel 4.7e-8), TREM2 2.1058310948785e-05 (5.6e-10), NECTIN2 0.000714815559569093 (9.3e-13); no NA
  component; the event log shows the two documented SNX22 Newton failures (Burden(1,25)-LINSIGHT: root 0.83 with
  residual 38.0; Burden(1,1)-LINSIGHT: root 20.6 with residual 1.53, both `newton_maxiter`) resolved by bisection, plus
  49 polished roots over the 3 sets.
* (b)(c) chunks 676 (SNX22), 326 (TREM2) and 849 (NECTIN2) are rerun with the unchanged driver logic into
  `revision/ukb/gw_coding_spa/rootcheck_validation/` and compared with the original chunk files by
  `scripts/compare_chunks.R` (rev 2), which attributes every logged event to its result row through the gene|category
  context and verifies the row's (#SNV, cMAC); disruptive_missense events on weight columns 0/13 are also attributed to
  the missense row's *-Disruptive columns. Per-column verdicts: a column without an event of its own must be identical
  to 1e-8 (check (c), inside event sets too — this is the real "normal branch unchanged" test, not the stored-p ≥ 0.05
  proxy of rev 1); a column touched only by Newton polishes may move by ≤ 1e-3 (b2); a column touched by a Newton
  failure / bisection / NA may change freely (b3). Gate criteria in addition: every event carries a context, maps to an
  existing row and column, matches the row's (#SNV, cMAC); every row has a `call_done` line and no call lacks a row;
  every event's p_normal < 0.05 (c2); no NA component; the three reference sets reproduce their values (e). Informative
  (not gating): the number of ambiguous (#SNV, cMAC) keys, call-order anomalies, disruptive→missense copies. Every
  changed set and component is listed with old/new values and verdicts in `results/validation/chunks_summary.txt`.
* `scripts/selftest_attribution.R` runs the first k genes of a chunk exactly as the driver does and lets
  `compare_chunks.R` compare them with the original chunk (gene filter): the end-to-end test of the attribution on real
  data without a full chunk rerun. **Run on 19 Oct for chunk 676, genes LACTB…SNX22 (12 genes, incl. APH1B whose
  disruptive call enters SPA, and SNX22): see section 11 and `results/validation_login/selftest_676/`.**
* The gate file is written only by `validate.R`; `jobs/03_release_gate.sbatch` (submitted with afterany on
  the validation) releases the scan only when it contains `"pass": true`, and every scan task re-checks it.

## 5. Genome-wide re-scan (jobs/03_gw_scan.sbatch)

`scripts/ukb_genome_wide_coding_scan.R` is a copy of `revision/scripts/ukb_gw_coding_spa_reconstructed.R` whose only
code change is the output directory (`gw_coding_spa/rootcheck/`, overridable through `UKB_ROOTCHECK_OUTDIR`); the library,
log names and resources are set in the sbatch. Outputs: `revision/ukb/gw_coding_spa/rootcheck/UKB_AD_Coding_<chunk>.Rdata`
and `rootcheck/rootlog/rootcheck_<chunk>.tsv` (solver events + one call_done line per burden call; a file with only
call_done lines means no solver event). Idempotent: an existing chunk output is skipped. Array 1-934 with throttle 60,
40 GB, 2 h, partition general.

Memory tiers (from the accounting of the original scan, job 4803341 and resubmissions 4804435/4806415/4843853): 909
chunks ran below 26 GB at 16 GB requested (median MaxRSS 9.7 GB, p90 17.5 GB); 25 OOM'd; chunks 13 20 70 281 283 329 345
405 needed 25–61 GB; chunks 142 282 816 OOM'd even at 64 GB and have no original result. The release script submits the
same sbatch three times: 1-934%60 at 40 GB (tasks of other tiers exit in seconds), the 8 heavy chunks at 64 GB and the 3
never-completed chunks at 180 GB/6 h (best effort). Any chunk still missing after the arrays is listed by the audit
(`results/audit/missing_chunks.txt`); resubmit it with
`sbatch --export=ALL,UKB_ROOTCHECK_TIER=<tier> --mem=<more> --array=<chunk> jobs/03_gw_scan.sbatch` (the tier must be the
chunk's tier: 40 unless it is in the lists above).

Cost: original 934 tasks, median 21 min, p90 38 min, max 76 min (≈ 390 task-hours; 25 OOM resubmissions on top). The
patched solver adds only the residual checks (a few O(n) passes per SPA component) and ~5 log appends per gene, so expect
the same ≈ 390 task-hours; with throttle 60, ≈ 7 h of wall time after the external dependencies clear. At 40 GB the
billing weight is ≈ 11 CPU-equivalents per task (TRESBillingWeights CPU=1, Mem=0.25G), ≈ 4,300 SU.

Post-scan audit (jobs/04_scan_audit.sbatch → results/audit/): `compare_chunks.R` over all chunks (every component and
aggregate p of ~87k sets, per-column verdicts), then `scan_audit_rootcheck.R`: count of sets whose STAAR-B changed by
more than 1e-6 relative, by category and magnitude; lists of sets crossing 1e-3, 1e-4, 1e-5 and 5e-7 in either
direction; the root-failure census from the event logs (components that fell back to bisection, polished roots,
bisection failures = NA in the output, by category, with the number of burden calls and of components that entered
SPA from the call_done lines as denominators); NA components; coverage; the same attribution/(b1)/(b2)/(b3)/(c)/(c2)
assertions as the validation; synonymous and functional calibration numbers old vs new (informative; the S9 builder is
the official one).

## 6. Downstream recomputation (jobs/05_downstream.sbatch, 4 tasks) and tables (jobs/06_tables.sbatch)

Copies of the original scripts, changed only in their output names (each copy's header lists the changes):

| copy (scripts/) | of (revision/scripts/) | writes |
|---|---|---|
| ukb_03_rootcheck_spa_lookups.R | ukb_03_reconstructed_spa_lookups.R | tables/review/ukb_rootcheck_genesets.csv, ukb_rootcheck_sv_lookup6.csv |
| ukb_single_variant_direct_test.R | ukb_sv_direct_reconstructed.R | tables/review/ukb_rootcheck_sv.csv (22 variants + APOE) |
| ukb_cond_rootcheck_spa.R | ukb_cond_reconstructed_spa.R | tables/review/ukb_conditional_rootcheck.csv |
| ukb_nectin2_driver_lookup_rootcheck.R | ukb_nectin2_driver_lookup.R | tables/review/ukb_nectin2_region_sv_rootcheck.csv |
| ukb_09_gw_calibration_rootcheck.R | ukb_09_gw_calibration_reconstructed.R | tables/supp_v3/S9_ukb_gw_reconstructed_SPA_burden_rootcheck.csv, S9b_…_hits_…_rootcheck.csv, S9c_…_focus_genes_…_rootcheck.csv, figures/supp_fig_S9_ukb_gw_riga_vs_noriga_SPA_burden_rootcheck.png |
| ukb_coding_gene_lookup.R | ukb_coding_lookup_reconstructed.R | tables/supp_v3/S10_ukb_coding_lookup_rootcheck.csv |
| ukb_build_replication_tables.py | ukb_build_S11_S13_reconstructed.py | tables/supp_v3/S11_ukb_geneset_replication_rootcheck.csv, S13_ukb_sv_direct_rootcheck.csv |
| compare_tables_rootcheck.py | — | results/tables/TABLE_COMPARISON_rootcheck.txt (+ _changes.csv; copied to tables/supp_v3/ukb_rootcheck_table_comparison.txt): every entry that differs from the current tables |

Note that the current S9 files live in `revision/tables/` (not supp_v3); the root-checked ones go to supp_v3 as asked.
The downstream solver event logs are in `results/downstream_rootlog/`.

## 7. Submission

```
bash SUBMIT.sh                 # build -> validate -> guarded release
bash SUBMIT.sh 5712345 5712346 # same, with the scan arrays queued afterany those jobs (e.g. the item-1 arrays)
```
The external ids travel to the release job COLON-separated (`--export=ALL,ROOTCHECK_EXTDEPS=id1:id2,ROOTCHECK_EXTDEPS_N=2`,
no whitespace inside an --export value); the release script splits on ':' and ABORTS (exit 4, nothing submitted) if the
value does not parse back to exactly ROOTCHECK_EXTDEPS_N numeric ids. Every sbatch in the release script is guarded
(`submit()` helper: a failed submission or an empty id aborts with exit 3 and is visible in the log/manifest). Job ids
are appended to `results/submission_manifest.txt`; the release job appends the ids of the arrays, downstream, audit and
tables. Logs: `${ADSP_ROOT}/log/review/rc_<name>_<jobid>[_<task>].out/.err`.
If the validation fails: nothing is released; read `VALIDATION.txt`, fix, re-run `sbatch jobs/02_validate.sbatch` and
then `bash jobs/03_release_gate.sbatch [external ids]` (it works as a plain shell script on the login node too).
Dry runs done on 19 Oct (rev 2): `sbatch --test-only` on every sbatch (all three scan tiers, and the release job with a
colon-separated --export value), `bash -n` on every shell script, R parse of every R file under `module load r/4.4.0`,
`py_compile` of both python files, and the release script's dependency parsing exercised without a gate file (2 ids
with N=2 → proceeds to the gate; N=3, a space-separated value, or a non-numeric argument → ABORT exit 4).

## 8. Files

```
README.md  SUBMIT.sh  build_patched_staar.sh  env_rootcheck.sh  VALIDATION.txt (written by the validation job)
src_untrusted/   the two upstream tarballs + extracted trees (unmodified; nothing is executed from here)
patches/         STAAR_saddlepoint_root_check.patch, STAARpipeline_saddlepoint_root_check.patch (rev 2)
build_work/      the patched trees the diffs were generated from (reference only)
build/           fresh patched trees + install logs + BUILD_RECORD.txt + BUILD_OK (written by build_patched_staar.sh)
jobs/            01_build, 02_validate, 03_release_gate, 03_gw_scan, 04_scan_audit, 05_downstream, 06_tables (.sbatch)
scripts/         ukb_genome_wide_coding_scan.R, unit_test.R, validate.R, compare_chunks.R,
                 ukb_common.R, scan_audit_rootcheck.R, selftest_attribution.R, the downstream/table copies,
                 compare_tables_rootcheck.py
gate/            ROOTCHECK_GATE.json (written by the validation job)
results/         validation/, validation_login/ (unit test, reference sets and the chunk-676 self-test run on 19 Oct),
                 audit/, downstream_rootlog/, tables/, submission_manifest.txt
```

## 9. Open questions / caveats

* Chunks 142, 282 and 816 never completed in the original scan (OOM at 64 GB; MaxRSS 67 GB when killed). The 180 GB
  tier is a guess; if they fail again the three chunks stay absent from both scans (the S9 calibration is computed on
  the chunks present, as before).
* The validation requires bit-identity for every column without an event of its own, inside event sets too. If any
  such column changes, the gate fails on purpose; the comparison lists it.
* `newton_polished` events are expected to be common (original Newton stops were accepted at |K1 − q| ≤ 1.2e-4): the
  affected p values move by ≤ ~1e-4 relative. The audit reports them separately from real root failures.
* The 2(a) component check against `components_guarded.csv` uses 1e-5 relative (two different root finders at 1e-6
  tolerance); the STAAR-B aggregates use 1e-6 as required.
* The release job runs `sbatch` from a compute node (interact). If the cluster refused that, run
  `bash jobs/03_release_gate.sbatch [ids]` on the login node after rc_validate has finished.
* The downstream copies keep the original scripts' behavior, including `ukb_09_gw_calibration_rootcheck.R` reading the
  noRIGA/RIGA null files for its figure annotation (as the original does) and the S13 builder's dependence on
  `supp_v3/S2_fdr_single_variants.csv`.
* Six R files of STAARpipeline now carry a Sys.setenv line before each burden call (diagnostic context). The gate
  proves them inert on real data (event-free columns bit-identical); the downstream scripts that use single-category
  functions get the context too, the noncoding functions (not used here) do not.
* p_normal in the log is the normal-approximation p of the component, and the original scan stored the SPA p for such
  components, so the two are not directly comparable; the reviewers' suggested cross-check is implemented as (c2):
  every event's p_normal < 0.05 (the component really entered SPA through the filter), plus the per-column identity test.
* Two behaviors of the ORIGINAL solver seen in the runs that will shape the genome-wide comparison: (i) when the
  original Newton hit max_iter with a Lugannani–Rice statistic below 38 it reset the root to 0 (p = 1) and the caller
  fell back to the OLD bisection, whose tolerance was 1.2e-4 absolute — such components (e.g. the real-null case H,
  the APH1B plof_ds/disruptive sets of chunk 676) change by ~1e-4 relative, i.e. they are "tiny changes" explained by
  newton_failure + bisection events, not by polishes; (ii) when the root lay outside the original ±100 bracket, the
  old bisection returned x = 0 → p = 1 and the component silently dropped out of the Cauchy combination; the patched
  bracket expansion finds those roots (APH1B plof_ds: roots at −275 … −3779, `bisection_ok_expanded`), so such
  components change from 1 to their real value — material changes, explained by the bisection events, expected
  genome-wide in sets with very rare high-weight variants.

## 10. Fixes applied after review (rev 2, 19 Oct)

1. **Disruptive call attribution** — events of a gene's disruptive_missense call are now attributed to its own row AND,
   for weight columns 0 (Burden(1,25)) and 13 (Burden(1,1)), to the missense row's `Burden(1,25)-Disruptive` /
   `Burden(1,1)-Disruptive` columns (`events_by_set()`, `attribution = "disruptive_copy"`); the missense aggregates are
   marked touched accordingly. (APH1B in chunk 676 is such a case; the self-test covers it.)
2. **n_weights** — the disruptive_missense call carries 26 weight columns like every other category (verified on chunk
   676: 29 p-columns = 26 + 3 aggregates); the 2-weight assumption came from the exact-mask login script (which calls
   the disruptive subset with annotation_phred = NULL, as the targeted run did) and is gone: n_weights is no longer part
   of any key, and the comment in `ukb_common.R` and README §2/§4 are corrected.
3. **Exact attribution** — the patched library records `<gene>|<category>` in every log line (STAARpipeline sets
   `STAAR_ROOTCHECK_CONTEXT` before each call; rootcheck.cpp copies it) and one `call_done` line per burden call. Events
   are attributed by context and verified against the row's (#SNV, cMAC); the number of ambiguous (#SNV, cMAC) keys is
   reported (`*_ambiguous_keys.csv`, status field `n_ambiguous_keys`) but no longer used. The call order is checked as
   a diagnostic (`n_call_order_anomalies`), with the correct code order plof_ds, plof, synonymous, missense,
   disruptive_missense.
4. **Check (c) per component** — within event sets every column without an event of its own must be identical to 1e-8
   (`verdict = VIOLATION_untouched_changed` otherwise); polish-only columns ≤ 1e-3; failure-touched columns free. The
   stored-p ≥ 0.05 proxy is gone. (c2) checks p_normal < 0.05 for every event.
5. **Unit test** — cases B4 (all four callers, incl. `STAAR_B_SPA_SMMAT_sparse`), G (no-root construction reaching every
   caller → NA in the component column) and H (the scan's path with the real sparse-kins null through
   `STAAR_Binary_SPA`, forced Newton failure and forced NA). The validation job passes the null file.
6. **SUBMIT.sh / release script** — colon-separated `ROOTCHECK_EXTDEPS` plus `ROOTCHECK_EXTDEPS_N`; the release script
   splits on ':', validates every token and aborts on a count mismatch.
7. **Release script guards** — every sbatch (scan tiers, downstream, audit, tables) goes through a guarded `submit()`
   helper; a failed submission or an empty id aborts with a non-zero exit before the manifest line is written.

## 11. Results of the rev-2 login-shell runs (19 Oct)

Unit test (`results/validation_login/unit_test_rev2.log`, re-run of the patched mode after a bookkeeping fix of one
check in `unit_test_rev2_patched_rerun.log`; CSVs `unit_test_checks_{original,patched}.csv`,
`unit_test_values_{original,patched}.csv`, event log `unit_test_rootlog.tsv`):
* ORIGINAL library (8/8 defect-reproducing checks): case A Newton root 0.554 with residual 14.4, p 4.8e-219 (reference
  0.03166); case B caller p 1.76e-5 vs reference 0.02831, and the SAME wrong value from `STAAR_B_SPA`,
  `STAAR_B_SPA_SMMAT` and `STAAR_B_SPA_SMMAT_sparse` (B4) and from STAARpipeline's single-variant copy (E); case G
  (no root reachable): every caller returns 1 (the component silently dropped); case H (real null): the original's
  partial safeguard happened to catch this construction (caller 0.0046044 / 0.0046035 vs reference 0.0046038 /
  0.0046037, 1.3e-4 / 4.3e-5 relative — the OLD bisection's looser tolerance), forced no-root → 1 everywhere.
* PATCHED library (25/25 checks after the fix; the first run failed only the check that expected 4 `call_done` lines
  where 5 burden calls had been made): A → NA, bisection root residual 9.8e-7, p 0.03166011 = reference to 1e-6; B →
  caller 0.0283110599 = bisection = reference; B4 → all four callers 0.0283110599 (1e-8), five `call_done` lines with
  n_spa 1, p_normal < 0.05 logged by the filtered callers; C and G → NA at the solver and in the component column of
  all four callers (log: newton_nonfinite 4, bisection_nobracket 4); D → verified root, no event; E → 0.0283110599;
  H → `STAAR_Binary_SPA` with the real null (`STAAR_B_SPA_SMMAT_sparse` path): Burden(1,25)-FAKE and Burden(1,1)-FAKE
  entered SPA, Newton NA (max_iter in both tails), bisection (two of them after bracket expansion), caller
  0.0046037812 / 0.0046037071 vs reference 0.0046037833 / 0.0046037058; `call_done` with n_variants 301, cMAC 302,
  n_weights 4, n_spa 2; forced no-root → NA, NA, NA, NA and STAAR-B 1, log bisection_nobracket ×4.

Reference sets with the rev-2 build (`unit_test_rev2.log`, section "refsets (rev-2 build)"; `refsets_rootlog.tsv` in the
17-column format): identical to the rev-1 values — SNX22 synonymous STAAR-B 0.0518593585025056 (rel 4.7e-8 to the
reference 0.0518593609546751), TREM2 missense 2.1058310948785e-05 (5.6e-10), NECTIN2 missense 0.000714815559569093
(9.3e-13); no NA component; max component difference to the guarded recomputation 4.0e-7 (SNX22), 7.5e-10, 3.4e-12. Event
log: 5 `call_done` lines (three full sets, two disruptive subsets), the two documented SNX22 Newton failures
(newton_maxiter at root 0.83/residual 38.0 and root 20.6/residual 1.53) plus two newton_unverified lower-tail roots, all
four resolved by bisection (one after bracket expansion), and 49 polished roots.

Attribution self-test on real data (`results/validation_login/selftest_676/`: chunk 676, genes LACTB RPS27L RAB8B APH1B
CA12 USP3 FBXL22 HERC1 DAPK2 CIAO2A SNX1 SNX22 run exactly as the driver, compared with the original chunk file by
`compare_chunks.R` with the gene filter; `cmp_summary.txt`, `cmp_status.json`, `cmp_events.csv`, `cmp_calls.csv`):
* Self-comparison first (`results/validation_login/selfcmp_676_partial/`, the reviewers' suggestion): the original
  chunk 676 compared with ITSELF, with the rerun's event log of the first four genes attached: 20 sets identical; 121
  events attributed by context (+4 disruptive→missense copies for APH1B, whose disruptive call enters SPA with
  Newton failures on weight columns 0 and 13), (#SNV, cMAC) verified for every event and every one of the 20
  `call_done` lines (76 SPA components, 0 calls without a row, 0 rows without a call, 0 call-order anomalies); 508
  untouched columns identical, 67 of them inside event sets; 80 touched columns; 1 ambiguous (#SNV, cMAC) key
  (LACTB/RPS27L disruptive_missense, 11/13) reported and unused; (b1)/(b2)/(b3)/(c)/(c2) all TRUE.
* Real rerun of the 12 genes vs the original chunk: RERUN_PLACEHOLDER
