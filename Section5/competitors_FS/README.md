# Section5/competitors_FS — competitor estimator for the Section 5 simulations

Contents
- `est_iv_fullsample.m`: full-sample IV (2SLS) estimator of beta, i.e. the rival benchmark
  \hat{beta}_{FS} reported in Tables 7 and 8 of the paper (no sub-sample selection).

Provenance
- Point estimate: literal port of `estimate_full_sample` / `OLS_Robust_full` of the original
  Section 5 replication package (same projection expression `Z*(Z'Z)^-1*Z'` and backslash
  solve, hence bit-exact with the published \hat{beta}_{FS}).
- Standard errors (not used in the paper's tables): White HC1 by default or Newey-West HAC
  (`Hac_option = true`, requires `m_break/NeweyWest.m` on the path), with structural residuals
  on the observed regressors, exactly as in `est_2sls_pilate.m`.
- With identical regressors in the two equations, full-sample SUR/FGLS coincides with
  equation-by-equation OLS (Zellner), so \hat{beta}_{FS} is also the full-sample counterpart of
  the FGLS estimator: one benchmark covers both.

Use
- Called only by `Section5/runSection5_Simul.m` (Stage 4) when its third argument
  `include_competitors` is `true` (default). It is not part of the standalone pi-LATE
  functions (`est_2sls_pilate.m`, `est_fgls_pilate.m`) nor of `Dynamic_LATE_Analysis.m`.
- `runSection5_Simul(n_sims, run_target, false)` reproduces exactly the previous outputs
  (same folder and file names, same `Paper_Results` fields, same tables).
- With `include_competitors = true` the output folder and file names carry the token `_FS`
  (`Paper_Exhibits_Results_<n>sims_SEC5_<TARGET>_FS`, `Detailed_Logs_..._FS_<timestamp>.txt`,
  `Paper_Exhibits_Suite_SEC5_..._FS_<timestamp>.mat`, `Paper_Exhibits_SEC5_<n>sims_FS_<timestamp>.txt`);
  `Paper_Results.<scenario>` gains the fields `Bias_FS` and `MSE_FS`; Tables 7 and 8 gain the
  rows Bias(\hat{beta}_{FS}) and MSE(\hat{beta}_{FS}) in the published order (FS, OLS, FGLS).
- DGP, seeds (`d_idx*1e6 + beta_idx*1e5 + i_sim`) and the pi-LATE stages are unchanged, so the
  2SLS and FGLS columns are identical with and without the competitor.

Verification
- Same seeds as the original package: the FS columns must coincide with the columns
  `Bias Full` / `MSE beta Full` of the original `Res_*.txt` result files (e.g. Tab. 6-7,
  rho = 0.25, pi0 = 0.6, i.i.d. errors: bias x100 = -11.995, -0.323, -0.162 and
  MSE x100 = 28.02, 1.176, 0.284 for d = 4, 16, 32).
