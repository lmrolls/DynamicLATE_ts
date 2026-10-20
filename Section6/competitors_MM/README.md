# Section6/competitors_MM — competitor tests for the Section 6 simulations

Contents (legacy code, unchanged)
- `computeMMtests2.m`, `computeMM_IVstats2.m`: Magnusson and Mavroeidis (2014) "MM" test
  statistics — GMM-S, split-S (GMM), qLL-S, exp-S, ave-S (first file) and the IV split-sample
  family split-S (IV), split-KLM, split-JKLM, split-CLR with the rank statistic rk-stat (second file).
- `getMagnussonAsymptoticCV.m`: simulated asymptotic critical values of qLL-S, exp-S, ave-S.
- `split_CLR_cv_seed.m`: simulated conditional critical value of split-CLR (per replication).
- `HACest.m`, `HACest_pw.m`, `HACest_tr.m`, `HAC_Sigma_hat.m`, `HACC_Sigma_hat.m`, `lag0.m`:
  HAC helpers from the original MM replication code.

Provenance
- All ten files are byte-identical (up to line endings) to the versions in the original
  Section 6 replication package (`Monte_Carlo/section_6` and `Magnusson_Mavroeidis/functions`).
- Known legacy conventions deliberately preserved in this branch: Newey-West kernel with lag
  floor(T^(1/3)); trimming 0.05 (the MM paper uses 0.15); the NW branch of split-KLM/JKLM/CLR,
  the form of the exp-S critical value and the scalar Omega of the homoskedastic branch are as in
  the legacy code (see CHANGELOG for the items to be revisited in the `fixes` branch).

Use
- Called only by `Section6/runSection6_Simul.m` when its third argument `include_competitors`
  is `true` (default). Not part of the standalone `robust_LATE_inference.m` nor of
  `Dynamic_LATE_Analysis.m`.
- `runSection6_Simul(n_sims, run_target, false)` reproduces exactly the previous pi-LATE-only
  outputs (same folder and file names, same `Paper_Results` fields with 6 columns, same tables).
- With `include_competitors = true` the output folder and file names carry the token `_MM`
  (`Paper_Exhibits_Results_<n>sims_SEC6_<TARGET>_MM`, `Detailed_Logs_..._MM_<timestamp>.txt`,
  `Paper_Exhibits_Suite_SEC6_..._MM_<timestamp>.mat`, `Paper_Exhibits_SEC6_<n>sims_MM_<timestamp>.txt`):
  `size` (Tables 9-10) and `power` (Figures 7-8) have 15 columns — AR_AMS, LM_AMS, LR_AMS, AR_ST,
  LM_ST, LR_ST, GMM_S_MM, split_S_MM_Feasible, qLL_S_MM, exp_S_MM, ave_S_MM, split_S_IV,
  split_KLM_IV, split_JKLM_IV, split_CLR_IV — plus `raw` with the means of M2_AMS, M2_ST, rk_stat_IV.
- Critical values: chi-square for GMM-S, split-S, split-KLM, split-JKLM; simulated asymptotic CVs
  for qLL-S, exp-S, ave-S (10000 draws, T = 10000); split-CLR CV simulated per replication with
  seed `rng_seed + 3` (legacy convention). For Table 9 the qLL-S/exp-S/ave-S critical values are
  re-simulated for every scenario after `rng(d*1000 + pi0*100 + Xactive)`, exactly as in the
  legacy script `simulFun_TAB_9.m` (switch `mm_cv_seed_per_scenario` inside the script; `false`
  keeps the once-per-run values). Figures 7-8 use size-adjusted power (empirical 95th percentile
  under H0) for all statistics.
- DGP, seeds and the pi-LATE stages are unchanged, so the six pi-LATE columns are identical with
  and without the competitors.

Verification
- Table 9 with `include_competitors = true` must reproduce the legacy Table 9 results
  (`results_TAB9_seed`), including the qLL-S, exp-S and ave-S columns thanks to the per-scenario
  seeding of their critical values; the six pi-LATE columns must coincide with the previous
  pi-LATE-only run.
