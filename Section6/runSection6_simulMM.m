function Paper_Results = runSection6_simulMM(n_sims_input, run_target)
% RUNSECTION6_SIMULMM Executes the Monte Carlo suite for Section 6 exhibits,
% EXTENDED to include the Magnusson & Mavroeidis (2014) "MM" test family
% (computeMMtests2.m, computeMM_IVstats2.m) alongside the existing AMS/ST
% (pi-LATE) tests.
%
% -------------------------------------------------------------------------
% DESIGN CONTRACT (read before modifying anything in this file):
%
%   1. This file is a STANDALONE COPY of runSection6_Simul.m, not a
%      refactor. runSection6_Simul.m is NEVER called, imported, or
%      modified by this file. Every line that reproduces AMS/ST logic
%      (generate_DGP, robust_LATE_inference calls, seed formulas, DGP
%      parameter grids) is a byte-for-byte copy of the corresponding line
%      in runSection6_Simul.m. Do not "clean up" or "optimize" those lines
%      (no inv()->\ changes, no restructuring) -- any such change breaks
%      the reproducibility guarantee this file exists to preserve.
%
%   2. MM statistics are ALWAYS computed here (no on/off flag). This file
%      is the "with MM" alternative to runSection6_Simul.m, not a
%      superset controlled by an option.
%
%   3. Seed discipline: every block that consumes randomness performs an
%      explicit rng(seed,"twister") at its own entry point and nowhere
%      else. This is already true of generate_DGP and of
%      LR_cv_legacy_exact (inside robust_LATE_inference, offsets +1/+2).
%      The MM addition reuses the SAME rng_seed already computed for
%      generate_DGP/robust_LATE_inference in that iteration, with two
%      consequences:
%        - computeMMtests2.m and computeMM_IVstats2.m draw NO randomness
%          internally (verified: no rand/randn/rng calls in either file),
%          so calling them does not disturb any stream.
%        - split_CLR_cv_seed.m (the one MM component that depends on
%          per-replication data, via rk_stat_IV) is seeded with
%          rng_seed+3, exactly continuing the legacy offset convention
%          already used by robust_LATE_inference for +1 (AMS) and +2 (ST).
%          Because every rng() call is a hard reset, this does not
%          interact with +1/+2 regardless of call order.
%      getMagnussonAsymptoticCV.m has no seed parameter and is called
%      exactly once per this file's execution (see point 5) -- this is a
%      known, accepted limitation, not a bug: its output is not tied to
%      any per-replication seed and is not claimed to be exactly
%      reproducible against an isolated legacy run of the same function.
%
%   4. Trimming: computeMMtests2.m/computeMM_IVstats2.m use trim=5%
%      (legacy convention), NOT the 15% specified in the Magnusson &
%      Mavroeidis (2014) supplemental appendix table captions. This is
%      LEFT UNCHANGED, per explicit instruction: fidelity to the legacy
%      codebase, not to the paper's stated tuning parameter.
%
%   5. Asymptotic/simulated critical values that do NOT depend on
%      per-replication data (GMM-S, split-S-Feasible, qLL-S, exp-S,
%      ave-S, split-S(IV), split-KLM, split-JKLM) are computed EXACTLY
%      ONCE for the whole run, before any loop, because
%      q_instr_MM=1, p_param_dim_MM=1, p_zeta_MM=0, alpha=0.05 are
%      constant across every scenario in Tables 9/10 in the legacy
%      scripts. Recomputing them per scenario would be functionally
%      redundant (10,000-draw simulations repeated for nothing) and is
%      NOT what the legacy scripts effectively needed, since none of
%      them varied these inputs.
%
%   6. split-CLR is the ONE MM statistic with a data-dependent critical
%      value (via rk_stat_IV per replication) and is therefore computed
%      per replication, inside the parfor, exactly mirroring how
%      LR_cv_legacy_exact works for CLR_T/CLR_T(S_hat).
%
%   7. rk_stat_IV has NO rejection rule anywhere in the legacy code
%      (simulFun_TAB_9.m explicitly excludes it from
%      all_test_names_stats and places it in
%      all_test_names_M2_components instead). It is reported as a raw
%      diagnostic value (mean, or mean trajectory over beta in
%      Figures 7/8), never as a rejection frequency. M2_AMS and M2_ST
%      receive the same treatment for consistency, since they play the
%      identical internal role (conditioning statistic for the CLR
%      critical value) and were already excluded from AMS/ST's own
%      rejection-frequency table in the legacy code.
%
%   8. In Figures 7/8, ALL statistics (AMS/ST and MM alike) use the
%      SAME size-adjusted-power mechanism already present in
%      runSection6_Simul.m: the empirical 95th percentile of the raw
%      statistic under H0 is used as the critical value (the
%      robust_LATE_inference cv_LR output is discarded there, exactly as
%      in the unmodified file). No asymptotic/simulated CV
%      (theoretical_CVs_MM) and no split_CLR_cv_seed call are used
%      inside the Figure 7/8 routines -- only inside Tables 9/10.
%
%   9. Naming: MATLAB field names for the MM statistics are the legacy
%      field names verbatim (GMM_S_MM, split_S_MM_Feasible, qLL_S_MM,
%      exp_S_MM, ave_S_MM, split_S_IV, split_KLM_IV, split_JKLM_IV,
%      split_CLR_IV, rk_stat_IV) -- these do not exist in the
%      Casini-McCloskey-Rolla-Pala paper and are not "translated" into
%      paper notation because the paper does not report them. Display
%      labels in the .txt export use the typographic convention of the
%      paper's own labels (hyphenated) applied consistently to the MM
%      names: GMM-S, split-S (GMM), split-S (IV), split-KLM,
%      split-JKLM, split-CLR, qLL-S, ave-S, exp-S, rk-stat.
%      AR_T is included in the export even though the paper's own
%      supplement states AR_T = LM_T for q=1 and therefore never reports
%      it separately -- included here because "report every statistic"
%      was an explicit instruction overriding paper-format filtering.
%
%  10. Both split-S variants (split_S_MM_Feasible from computeMMtests2.m
%      and split_S_IV from computeMM_IVstats2.m) are reported as SEPARATE
%      columns. The paper has a single "split-S" column; which of the two
%      legacy families it corresponds to was never resolved (the two
%      helpers implement genuinely different constructions -- GMM-moment
%      based vs. IV/Kleibergen based). Reporting both avoids an
%      unjustified choice between them.
%
%  11. KNOWN INHERITED DEFECT, deliberately not fixed here: the exp-S
%      critical value produced by getMagnussonAsymptoticCV.m ("Form 1")
%      does not use the same scaling of the bridge functional as the
%      exp_S_MM statistic computed in computeMMtests2.m. This was
%      verified analytically (Jensen's inequality on the FCLT limit of
%      the split-sample statistic) and confirmed by direct Monte Carlo
%      (empirical size ~0.032 vs nominal 0.05 under a stylized DGP).
%      Left untouched per explicit instruction: legacy fidelity, no
%      corrections.
%
%  12. params_MM_base / params_MM_IV_base (the static field set: hac
%      method, trim proportion, instrument count, Newey-West lag, field
%      name lists) are built EXACTLY ONCE per routine block (Table 9,
%      Table 10, Figure 7, Figure 8), immediately after T_val is fixed
%      for that routine -- not once per parfor iteration. Inside each
%      parfor iteration, only Mx_s_local_full_T is added, via
%      params_MM_it = params_MM_base; params_MM_it.Mx_s_local_full_T = Mx_MM;
%      This mirrors the legacy pattern (base struct copied, one field
%      overwritten) rather than rebuilding every field from scratch on
%      every replication.
% -------------------------------------------------------------------------

    clc;
    restoredefaultpath;
    addpath('..\m_break');      % Newey-West kernels (NeweyWest, NeweyWest_N1N2) -- unchanged path, verbatim
    addpath('helpersMM');       % computeMMtests2, computeMM_IVstats2, getMagnussonAsymptoticCV, split_CLR_cv_seed
    warning('off','all');

    % --- OPTIONAL INPUTS AND DEFAULT HANDLING (verbatim) ---
    if nargin < 1 || isempty(n_sims_input), n_sims = 10000; else, n_sims = n_sims_input; end
    if nargin < 2 || isempty(run_target), run_target = 'all'; end

    valid_targets = {'all', 'tab9', 'tab10', 'fig7', 'fig8'};
    if ~ismember(lower(run_target), valid_targets)
        error('Error: run_target must be ''all'', ''tab9'', ''tab10'', ''fig7'', or ''fig8''.');
    end

    run_tab9  = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab9');
    run_tab10 = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab10');
    run_fig7  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig7');
    run_fig8  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig8');

    % --- DESCRIPTIVE STRINGS FOR LOGS AND DIRECTORY SETUP (MM token added) ---
    folder_name = sprintf('Paper_Exhibits_Results_%dsims_SEC6_%s_MM', n_sims, upper(run_target));
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end

    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    timestamp_disp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    log_filename = fullfile(folder_name, sprintf('Detailed_Logs_%dsims_SEC6_%s_MM_%s.txt', n_sims, upper(run_target), timestamp_file));
    fid_log = fopen(log_filename, 'w');
    if fid_log == -1, error('Unable to create the log file.'); end

    fprintf_both(fid_log, '===================================================================\n');
    fprintf_both(fid_log, ' STARTING SIMULATION SUITE: SECTION 6 (DYNAMIC LATE) -- MM-EXTENDED\n');
    fprintf_both(fid_log, ' Execution Target: %s\n', upper(run_target));
    fprintf_both(fid_log, ' Date/Time       : %s\n', timestamp_disp);
    fprintf_both(fid_log, ' Iterations      : %d\n', n_sims);
    fprintf_both(fid_log, ' Focus           : AMS, ST (pi-LATE) AND Magnusson-Mavroeidis (MM) Tests\n');
    fprintf_both(fid_log, ' Engine          : robust_LATE_inference (unmodified) + computeMMtests2 / computeMM_IVstats2 (helpersMM)\n');
    fprintf_both(fid_log, '===================================================================\n\n');

    Paper_Results = struct();

    % Core logic parameters (legacy math framework structurally embedded) -- VERBATIM, identical to runSection6_Simul.m
    params_core = struct('q', 1, 'max_breaks', 5, 'eps1', 0.05, 'vector_piL', (0.6:0.1:1), 'step_pi', 0.05, 'Hac_option', true);

    % --- Canonical ordered statistic names, used only for reporting/printing ---
    rej_names_all = {'AR_AMS','LM_AMS','LR_AMS','AR_ST','LM_ST','LR_ST', ...
                      'GMM_S_MM','split_S_MM_Feasible','qLL_S_MM','exp_S_MM','ave_S_MM', ...
                      'split_S_IV','split_KLM_IV','split_JKLM_IV','split_CLR_IV'};   % 15 columns, fixed order
    raw_names_all = {'M2_AMS','M2_ST','rk_stat_IV'};                                 % 3 columns, fixed order

    % =========================================================================
    % PRECOMPUTE MM ASYMPTOTIC / SIMULATED CRITICAL VALUES ONCE (see contract #5)
    % =========================================================================
    alpha_test_MM   = 0.05;
    q_instr_MM      = 1;
    p_param_dim_MM  = 1;
    p_zeta_MM       = 0;
    N_cv_sims_MM    = 10000;
    T_cv_approx_MM  = 10000;
    trim_prop_MM    = 0.05;    % legacy value, NOT the paper's 15% -- see contract #4

    fprintf_both(fid_log, '\n>>> Precomputing MM asymptotic/simulated critical values (once, alpha=%.2f)...\n', alpha_test_MM);

    theoretical_CVs_MM = struct();
    theoretical_CVs_MM.GMM_S_MM = chi2inv(1 - alpha_test_MM, max(1, q_instr_MM - p_zeta_MM));

    df_split_s_mm = 2*q_instr_MM - p_zeta_MM;
    if df_split_s_mm <= 0
        theoretical_CVs_MM.split_S_MM_Feasible = inf;
    else
        theoretical_CVs_MM.split_S_MM_Feasible = chi2inv(1 - alpha_test_MM, df_split_s_mm);
    end

    theoretical_CVs_MM.qLL_S_MM = getMagnussonAsymptoticCV('qLL-S', q_instr_MM, alpha_test_MM, N_cv_sims_MM, T_cv_approx_MM, trim_prop_MM, p_zeta_MM);
    theoretical_CVs_MM.exp_S_MM = getMagnussonAsymptoticCV('exp-S', q_instr_MM, alpha_test_MM, N_cv_sims_MM, T_cv_approx_MM, trim_prop_MM, p_zeta_MM);
    theoretical_CVs_MM.ave_S_MM = getMagnussonAsymptoticCV('ave-S', q_instr_MM, alpha_test_MM, N_cv_sims_MM, T_cv_approx_MM, trim_prop_MM, p_zeta_MM);

    theoretical_CVs_MM.split_S_IV   = chi2inv(1 - alpha_test_MM, 2*q_instr_MM);
    theoretical_CVs_MM.split_KLM_IV = chi2inv(1 - alpha_test_MM, p_param_dim_MM);

    df_jklm = 2*q_instr_MM - p_param_dim_MM;
    if df_jklm > 0
        theoretical_CVs_MM.split_JKLM_IV = chi2inv(1 - alpha_test_MM, df_jklm);
    else
        theoretical_CVs_MM.split_JKLM_IV = inf;
    end
    % split_CLR_IV has NO entry here: it is data-dependent, computed per replication via split_CLR_cv_seed (contract #6)

    fprintf_both(fid_log, '    GMM-S=%.4f  split-S(GMM)=%.4f  qLL-S=%.4f  exp-S=%.4f  ave-S=%.4f\n', ...
        theoretical_CVs_MM.GMM_S_MM, theoretical_CVs_MM.split_S_MM_Feasible, theoretical_CVs_MM.qLL_S_MM, theoretical_CVs_MM.exp_S_MM, theoretical_CVs_MM.ave_S_MM);
    fprintf_both(fid_log, '    split-S(IV)=%.4f  split-KLM=%.4f  split-JKLM=%.4f\n\n', ...
        theoretical_CVs_MM.split_S_IV, theoretical_CVs_MM.split_KLM_IV, theoretical_CVs_MM.split_JKLM_IV);

    % =========================================================================
    % ROUTINE: TABLE 9 (Empirical Size, i.i.d., T=400) -- AMS/ST block verbatim, MM appended
    % =========================================================================
    if run_tab9
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLE 9 (Empirical Size, i.i.d., T=400) -- MM-EXTENDED\n');
        fprintf_both(fid_log, '===================================================================\n');

        T_val = 400; beta_H0 = 0.0; rho_eu = 0.25;
        d_grid = [4, 8, 12, 16]; pi0_grid = [0.4, 0.6, 0.8]; Xactive_grid = [0, 1];

        % --- params_MM_base / params_MM_IV_base built ONCE for this routine (contract #12) ---
        params_MM_base.hac_method_MM = 'NW';
        params_MM_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_base.q_instr_MM = q_instr_MM;
        params_MM_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_base.test_names_MM_fields = {'GMM_S_MM','split_S_MM_Feasible','qLL_S_MM','exp_S_MM','ave_S_MM'};

        params_MM_IV_base.p_param_dim = p_param_dim_MM;
        params_MM_IV_base.q_instr = q_instr_MM;
        params_MM_IV_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_IV_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_IV_base.hac_method_MM = 'NW';

        for pi0_val = pi0_grid
            for Xactive = Xactive_grid
                for d_val = d_grid
                    run_label = sprintf('Tab9_d%d_pi0%02d_X%d', d_val, pi0_val*100, Xactive);
                    fprintf_both(fid_log, '\n--- Starting Iteration: Tab 9 [d=%d, pi0=%.1f, Xactive=%d] (MM-extended) ---\n', d_val, pi0_val, Xactive);

                    theta_H0 = d_val / sqrt(T_val);
                    rej_matrix = zeros(n_sims, 15);
                    raw_matrix = zeros(n_sims, 3);

                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = d_val*1e8 + (pi0_val*100)*1e6 + Xactive*1e5 + i_sims;   % IDENTICAL to runSection6_Simul.m

                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H0, 0, beta_H0, 1, rho_eu, Xactive, rng_seed);       % IDENTICAL call
                        [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);      % IDENTICAL call

                        rej_row_1_6 = [ ...
                            (stats_AMS.AR * params_core.q) > cv_LR.cv_AR, stats_AMS.LM > cv_LR.cv_LM, stats_AMS.LR > cv_LR.cv_AMS, ...
                            (stats_ST.AR * params_core.q)  > cv_LR.cv_AR, stats_ST.LM  > cv_LR.cv_LM, stats_ST.LR  > cv_LR.cv_ST ...
                        ];                                                                                                                 % IDENTICAL 6-vector
                        raw_row_1_2 = [stats_AMS.M2, stats_ST.M2];

                        % --- MM ADDITION (new; does not alter anything above) ---
                        Mx_MM = eye(T_val);
                        if Xactive == 1 && rank(X_raw'*X_raw) == size(X_raw,2)
                            Mx_MM = eye(T_val) - X_raw*((X_raw'*X_raw)\X_raw');
                        end

                        params_MM_it = params_MM_base;
                        params_MM_it.Mx_s_local_full_T = Mx_MM;
                        params_MM_IV_it = params_MM_IV_base;
                        params_MM_IV_it.Mx_s_local_full_T = Mx_MM;

                        res_MM    = computeMMtests2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_it);
                        res_MM_IV = computeMM_IVstats2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_IV_it);

                        rej_row_7_11 = [ ...
                            res_MM.GMM_S_MM > theoretical_CVs_MM.GMM_S_MM, ...
                            res_MM.split_S_MM_Feasible > theoretical_CVs_MM.split_S_MM_Feasible, ...
                            res_MM.qLL_S_MM > theoretical_CVs_MM.qLL_S_MM, ...
                            res_MM.exp_S_MM > theoretical_CVs_MM.exp_S_MM, ...
                            res_MM.ave_S_MM > theoretical_CVs_MM.ave_S_MM ...
                        ];

                        cv_splitCLR_it = NaN;
                        if ~isnan(res_MM_IV.rk_stat_IV)
                            cv_splitCLR_it = split_CLR_cv_seed(res_MM_IV.rk_stat_IV, p_param_dim_MM, q_instr_MM, alpha_test_MM, rng_seed+3);
                        end
                        rej_row_12_15 = [ ...
                            res_MM_IV.split_S_IV > theoretical_CVs_MM.split_S_IV, ...
                            res_MM_IV.split_KLM_IV > theoretical_CVs_MM.split_KLM_IV, ...
                            res_MM_IV.split_JKLM_IV > theoretical_CVs_MM.split_JKLM_IV, ...
                            res_MM_IV.split_CLR_IV > cv_splitCLR_it ...
                        ];
                        raw_row_3 = res_MM_IV.rk_stat_IV;

                        % Single sliced write per matrix, uniform (i_sims, :) subscript -- required by parfor classification
                        rej_matrix(i_sims, :) = [rej_row_1_6, rej_row_7_11, rej_row_12_15];
                        raw_matrix(i_sims, :) = [raw_row_1_2, raw_row_3];
                    end

                    emp_size  = mean(rej_matrix, 1, 'omitnan');
                    raw_means = mean(raw_matrix, 1, 'omitnan');
                    Paper_Results.(run_label) = struct('params', [d_val, pi0_val, Xactive], 'size', emp_size, 'raw', raw_means);

                    print_partial_results_MM(fid_log, emp_size, raw_means, rej_names_all, raw_names_all);
                end
            end
        end
    end

    % =========================================================================
    % ROUTINE: TABLE 10 (Empirical Size, i.i.d., T=200) -- AMS/ST block verbatim, MM appended
    % =========================================================================
    if run_tab10
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLE 10 (Empirical Size, i.i.d. errors, T=200) -- MM-EXTENDED\n');
        fprintf_both(fid_log, '>>> FORCING COMMON RANDOM NUMBERS (LEGACY SEED MATCH)\n');
        fprintf_both(fid_log, '===================================================================\n');

        T_val = 200; beta_H0 = 0.0; pi0_val = 0.6;
        rho_grid = [0.25, 0.50, 0.75]; d_grid = [4, 10, 16]; Xactive_grid = [0, 1];

        % --- params_MM_base / params_MM_IV_base built ONCE for this routine (contract #12) ---
        params_MM_base.hac_method_MM = 'NW';
        params_MM_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_base.q_instr_MM = q_instr_MM;
        params_MM_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_base.test_names_MM_fields = {'GMM_S_MM','split_S_MM_Feasible','qLL_S_MM','exp_S_MM','ave_S_MM'};

        params_MM_IV_base.p_param_dim = p_param_dim_MM;
        params_MM_IV_base.q_instr = q_instr_MM;
        params_MM_IV_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_IV_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_IV_base.hac_method_MM = 'NW';

        for Xactive = Xactive_grid
            for rho_val = rho_grid
                for d_val = d_grid
                    run_label = sprintf('Tab10_X%d_rho%02d_d%d', Xactive, rho_val*100, d_val);
                    fprintf_both(fid_log, '\n--- Starting Iteration: Tab 10 [Xactive=%d, rho=%.2f, d=%d] (MM-extended) ---\n', Xactive, rho_val, d_val);

                    theta_H0 = d_val / sqrt(T_val);
                    rej_matrix = zeros(n_sims, 15);
                    raw_matrix = zeros(n_sims, 3);

                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = 1*1e6 + i_sims*100 + 16e7;                                                                              % IDENTICAL to runSection6_Simul.m

                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H0, 0, beta_H0, 1, rho_val, Xactive, rng_seed);      % IDENTICAL call
                        [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);      % IDENTICAL call

                        rej_row_1_6 = [ ...
                            (stats_AMS.AR * params_core.q) > cv_LR.cv_AR, stats_AMS.LM > cv_LR.cv_LM, stats_AMS.LR > cv_LR.cv_AMS, ...
                            (stats_ST.AR * params_core.q)  > cv_LR.cv_AR, stats_ST.LM  > cv_LR.cv_LM, stats_ST.LR  > cv_LR.cv_ST ...
                        ];                                                                                                                 % IDENTICAL 6-vector
                        raw_row_1_2 = [stats_AMS.M2, stats_ST.M2];

                        % --- MM ADDITION ---
                        Mx_MM = eye(T_val);
                        if Xactive == 1 && rank(X_raw'*X_raw) == size(X_raw,2)
                            Mx_MM = eye(T_val) - X_raw*((X_raw'*X_raw)\X_raw');
                        end

                        params_MM_it = params_MM_base;
                        params_MM_it.Mx_s_local_full_T = Mx_MM;
                        params_MM_IV_it = params_MM_IV_base;
                        params_MM_IV_it.Mx_s_local_full_T = Mx_MM;

                        res_MM    = computeMMtests2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_it);
                        res_MM_IV = computeMM_IVstats2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_IV_it);

                        rej_row_7_11 = [ ...
                            res_MM.GMM_S_MM > theoretical_CVs_MM.GMM_S_MM, ...
                            res_MM.split_S_MM_Feasible > theoretical_CVs_MM.split_S_MM_Feasible, ...
                            res_MM.qLL_S_MM > theoretical_CVs_MM.qLL_S_MM, ...
                            res_MM.exp_S_MM > theoretical_CVs_MM.exp_S_MM, ...
                            res_MM.ave_S_MM > theoretical_CVs_MM.ave_S_MM ...
                        ];

                        cv_splitCLR_it = NaN;
                        if ~isnan(res_MM_IV.rk_stat_IV)
                            cv_splitCLR_it = split_CLR_cv_seed(res_MM_IV.rk_stat_IV, p_param_dim_MM, q_instr_MM, alpha_test_MM, rng_seed+3);
                        end
                        rej_row_12_15 = [ ...
                            res_MM_IV.split_S_IV > theoretical_CVs_MM.split_S_IV, ...
                            res_MM_IV.split_KLM_IV > theoretical_CVs_MM.split_KLM_IV, ...
                            res_MM_IV.split_JKLM_IV > theoretical_CVs_MM.split_JKLM_IV, ...
                            res_MM_IV.split_CLR_IV > cv_splitCLR_it ...
                        ];
                        raw_row_3 = res_MM_IV.rk_stat_IV;

                        % Single sliced write per matrix, uniform (i_sims, :) subscript -- required by parfor classification
                        rej_matrix(i_sims, :) = [rej_row_1_6, rej_row_7_11, rej_row_12_15];
                        raw_matrix(i_sims, :) = [raw_row_1_2, raw_row_3];
                    end

                    emp_size  = mean(rej_matrix, 1, 'omitnan');
                    raw_means = mean(raw_matrix, 1, 'omitnan');
                    Paper_Results.(run_label) = struct('params', [Xactive, rho_val, d_val], 'size', emp_size, 'raw', raw_means);

                    print_partial_results_MM(fid_log, emp_size, raw_means, rej_names_all, raw_names_all);
                end
            end
        end
    end

    % =========================================================================
    % ROUTINE: FIGURE 7 (Size-Adjusted Power Curve, i.i.d., T=200) -- AMS/ST block verbatim, MM appended
    % =========================================================================
    if run_fig7
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: FIGURE 7 (Size-Adjusted Power, T=200) -- MM-EXTENDED\n');
        fprintf_both(fid_log, '===================================================================\n');

        T_val = 200; beta_H0 = 0.0; pi0_val = 0.6; Xactive = 0;
        rho_grid = [0.25, 0.75]; d_grid = [4, 10, 16]; beta_grid = 0:0.1:1;

        % --- params_MM_base / params_MM_IV_base built ONCE for this routine (contract #12) ---
        params_MM_base.hac_method_MM = 'NW';
        params_MM_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_base.q_instr_MM = q_instr_MM;
        params_MM_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_base.test_names_MM_fields = {'GMM_S_MM','split_S_MM_Feasible','qLL_S_MM','exp_S_MM','ave_S_MM'};

        params_MM_IV_base.p_param_dim = p_param_dim_MM;
        params_MM_IV_base.q_instr = q_instr_MM;
        params_MM_IV_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_IV_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_IV_base.hac_method_MM = 'NW';

        for rho_val = rho_grid
            for d_val = d_grid
                run_label = sprintf('Fig7_rho%02d_d%d', rho_val*100, d_val);
                fprintf_both(fid_log, '\n--- Starting Power Curve: Fig 7 [rho=%.2f, d=%d] (MM-extended) ---\n', rho_val, d_val);

                theta_H1_strong = d_val / sqrt(T_val);
                stats_H0_matrix = zeros(n_sims, 15);

                fprintf_both(fid_log, '    > Computing empirical CVs under H0 (all 15 statistics, size-adjusted)...\n');
                parfor i_sims = 1:n_sims
                    warning('off','all');
                    rng_seed = 1*1e6 + i_sims*100 + 16e7;                                                                                  % IDENTICAL

                    [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, 0, beta_H0, 1, rho_val, Xactive, rng_seed);   % IDENTICAL
                    [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);              % IDENTICAL (cv_LR discarded, as in the unmodified file)

                    Mx_MM = eye(T_val);   % Xactive==0 throughout Fig 7 -> Mx_MM stays eye(T_val); branch kept for generality/robustness of this new code only
                    if Xactive == 1 && rank(X_raw'*X_raw) == size(X_raw,2)
                        Mx_MM = eye(T_val) - X_raw*((X_raw'*X_raw)\X_raw');
                    end

                    params_MM_it = params_MM_base;
                    params_MM_it.Mx_s_local_full_T = Mx_MM;
                    params_MM_IV_it = params_MM_IV_base;
                    params_MM_IV_it.Mx_s_local_full_T = Mx_MM;

                    res_MM    = computeMMtests2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_it);
                    res_MM_IV = computeMM_IVstats2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_IV_it);

                    stats_H0_matrix(i_sims, :) = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR, ...
                        res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                        res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                end
                emp_CVs = prctile(stats_H0_matrix, 95, 1);   % 15-column size-adjusted CV, same mechanism as AMS/ST (see contract #8)

                power_matrix = zeros(length(beta_grid), 15);
                raw_matrix_fig = zeros(length(beta_grid), 3);
                fprintf_both(fid_log, '    > Computing Power over beta grid (15 rejection-based + 3 raw-trajectory statistics)...\n');
                for b_idx = 1:length(beta_grid)
                    beta_alt = beta_grid(b_idx);
                    rej_H1 = zeros(n_sims, 15);
                    raw_H1 = zeros(n_sims, 3);
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = b_idx*1e7 + i_sims*100 + 25e7;                                                                            % IDENTICAL

                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, 0, beta_alt, 1, rho_val, Xactive, rng_seed); % IDENTICAL
                        [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);             % IDENTICAL

                        Mx_MM = eye(T_val);
                        if Xactive == 1 && rank(X_raw'*X_raw) == size(X_raw,2)
                            Mx_MM = eye(T_val) - X_raw*((X_raw'*X_raw)\X_raw');
                        end

                        params_MM_it = params_MM_base;
                        params_MM_it.Mx_s_local_full_T = Mx_MM;
                        params_MM_IV_it = params_MM_IV_base;
                        params_MM_IV_it.Mx_s_local_full_T = Mx_MM;

                        res_MM    = computeMMtests2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_it);
                        res_MM_IV = computeMM_IVstats2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_IV_it);

                        vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR, ...
                            res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                            res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                        rej_H1(i_sims, :) = vals > emp_CVs;
                        raw_H1(i_sims, :) = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                    end
                    power_matrix(b_idx, :) = mean(rej_H1, 1, 'omitnan');
                    raw_matrix_fig(b_idx, :) = mean(raw_H1, 1, 'omitnan');
                end
                Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix, 'raw', raw_matrix_fig);
                fprintf_both(fid_log, '    > Power computed successfully.\n');
            end
        end
    end

    % =========================================================================
    % ROUTINE: FIGURE 8 (Perverse Identification, i.i.d., T=200) -- AMS/ST block verbatim, MM appended
    % =========================================================================
    if run_fig8
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: FIGURE 8 (Perverse Identification, T=200) -- MM-EXTENDED\n');
        fprintf_both(fid_log, '===================================================================\n');

        T_val = 200; beta_H0 = 0.0; rho_val = 0.25; Xactive = 0;
        pi0_grid = [0.6, 0.8]; d_grid = [24, 32]; beta_grid = 0:0.1:1;
        perverse_types = {'local', 'fixed'};

        % --- params_MM_base / params_MM_IV_base built ONCE for this routine (contract #12) ---
        params_MM_base.hac_method_MM = 'NW';
        params_MM_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_base.q_instr_MM = q_instr_MM;
        params_MM_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_base.test_names_MM_fields = {'GMM_S_MM','split_S_MM_Feasible','qLL_S_MM','exp_S_MM','ave_S_MM'};

        params_MM_IV_base.p_param_dim = p_param_dim_MM;
        params_MM_IV_base.q_instr = q_instr_MM;
        params_MM_IV_base.trim_prop_cv_mm = trim_prop_MM;
        params_MM_IV_base.Lag_NW_s_local = floor(T_val^(1/3));
        params_MM_IV_base.hac_method_MM = 'NW';

        for p_idx = 1:length(perverse_types)
            p_type = perverse_types{p_idx};
            for pi0_val = pi0_grid
                for d_val = d_grid

                    theta_H1_strong = d_val / sqrt(T_val);
                    if strcmp(p_type, 'local'), theta_H1_weak = -0.5 * d_val / sqrt(T_val);
                    else, theta_H1_weak = -0.5; end

                    run_label = sprintf('Fig8_%s_pi0%02d_d%d', p_type, pi0_val*100, d_val);
                    fprintf_both(fid_log, '\n--- Starting Power Curve: Fig 8 [Type=%s, pi0=%.1f, d=%d] (MM-extended) ---\n', p_type, pi0_val, d_val);

                    stats_H0_matrix = zeros(n_sims, 15);
                    fprintf_both(fid_log, '    > Computing empirical CVs under H0 (all 15 statistics, size-adjusted)...\n');
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = 1*1e6 + i_sims*100 + 16e7;                                                                                % IDENTICAL

                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, theta_H1_weak, beta_H0, 1, rho_val, Xactive, rng_seed); % IDENTICAL
                        [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);               % IDENTICAL

                        Mx_MM = eye(T_val);
                        if Xactive == 1 && rank(X_raw'*X_raw) == size(X_raw,2)
                            Mx_MM = eye(T_val) - X_raw*((X_raw'*X_raw)\X_raw');
                        end

                        params_MM_it = params_MM_base;
                        params_MM_it.Mx_s_local_full_T = Mx_MM;
                        params_MM_IV_it = params_MM_IV_base;
                        params_MM_IV_it.Mx_s_local_full_T = Mx_MM;

                        res_MM    = computeMMtests2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_it);
                        res_MM_IV = computeMM_IVstats2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_IV_it);

                        stats_H0_matrix(i_sims, :) = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR, ...
                            res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                            res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                    end
                    emp_CVs = prctile(stats_H0_matrix, 95, 1);

                    power_matrix = zeros(length(beta_grid), 15);
                    raw_matrix_fig = zeros(length(beta_grid), 3);
                    fprintf_both(fid_log, '    > Computing Power over beta grid (15 rejection-based + 3 raw-trajectory statistics)...\n');
                    for b_idx = 1:length(beta_grid)
                        beta_alt = beta_grid(b_idx);
                        rej_H1 = zeros(n_sims, 15);
                        raw_H1 = zeros(n_sims, 3);
                        parfor i_sims = 1:n_sims
                            warning('off','all');
                            rng_seed = b_idx*1e7 + i_sims*100 + 25e7;                                                                          % IDENTICAL

                            [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, theta_H1_weak, beta_alt, 1, rho_val, Xactive, rng_seed); % IDENTICAL
                            [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);           % IDENTICAL

                            Mx_MM = eye(T_val);
                            if Xactive == 1 && rank(X_raw'*X_raw) == size(X_raw,2)
                                Mx_MM = eye(T_val) - X_raw*((X_raw'*X_raw)\X_raw');
                            end

                            params_MM_it = params_MM_base;
                            params_MM_it.Mx_s_local_full_T = Mx_MM;
                            params_MM_IV_it = params_MM_IV_base;
                            params_MM_IV_it.Mx_s_local_full_T = Mx_MM;

                            res_MM    = computeMMtests2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_it);
                            res_MM_IV = computeMM_IVstats2(Y_s, D_s, Z_raw, X_raw, T_val, beta_H0, params_MM_IV_it);

                            vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR, ...
                                res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                                res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                            rej_H1(i_sims, :) = vals > emp_CVs;
                            raw_H1(i_sims, :) = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                        end
                        power_matrix(b_idx, :) = mean(rej_H1, 1, 'omitnan');
                        raw_matrix_fig(b_idx, :) = mean(raw_H1, 1, 'omitnan');
                    end
                    Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix, 'raw', raw_matrix_fig);
                    fprintf_both(fid_log, '    > Power computed successfully.\n');
                end
            end
        end
    end

    % --- FINAL SHUTDOWN AND EXPORT ---
    reset_parallel_pool(fid_log);
    fclose(fid_log);

    out_filename_mat = fullfile(folder_name, sprintf('Paper_Exhibits_Suite_SEC6_%dsims_%s_MM_%s.mat', n_sims, upper(run_target), timestamp_file));
    save(out_filename_mat, 'Paper_Results');

    export_paper_exhibits_sec6_MM_to_txt(Paper_Results, n_sims, folder_name, run_target, timestamp_file, timestamp_disp, rej_names_all, raw_names_all);

    fprintf('\nExecution Completed Successfully (Target: %s, MM-extended). Outputs saved in: %s\\\n', upper(run_target), folder_name);
end


%% ========================================================================
%  LOCAL UTILITIES (Data Generation) -- verbatim copies from runSection6_Simul.m
%  MATLAB does not allow calling another function file's local functions
%  from outside that file, so these MUST be duplicated here unchanged.
% =========================================================================

function reset_parallel_pool(fid)
    fprintf_both(fid, '\n[SYSTEM] Executing memory cleanup and resetting Parallel Pool...\n');
    poolobj = gcp('nocreate');
    if ~isempty(poolobj), delete(poolobj); end
    num_workers = max(1, feature('numcores') - 2);
    parpool('local', num_workers);
    fprintf_both(fid, '[SYSTEM] Parallel Pool started with %d workers.\n', num_workers);
end

function fprintf_both(fid, fmt, varargin)
    fprintf(1, fmt, varargin{:});
    if fid ~= -1, fprintf(fid, fmt, varargin{:}); end
end

function print_partial_results_MM(fid, emp_size, raw_means, rej_names, raw_names)
    fprintf_both(fid, '    > Partial Rejection Rates (AMS/ST + MM, 15 statistics):\n');
    fprintf_both(fid, '      %s\n', strjoin(rej_names, ' | '));
    fprintf_both(fid, '      %s\n', strjoin(arrayfun(@(x) sprintf('%.4f',x), emp_size, 'UniformOutput', false), ' | '));
    fprintf_both(fid, '    > Raw diagnostic means (no rejection rule, %s):\n', strjoin(raw_names, ', '));
    fprintf_both(fid, '      %s\n\n', strjoin(arrayfun(@(x) sprintf('%.4f',x), raw_means, 'UniformOutput', false), ' | '));
end

function [Y, D, Z, X] = generate_DGP(T, pi0, theta_strong, theta_weak, beta, gamma, rho, Xactive, rng_seed)
    % VERBATIM copy of generate_DGP from runSection6_Simul.m -- do not alter.
    rng(rng_seed, "twister");
    Sigma_eu = [1 rho; rho 1];
    Z = 1 + randn(T, 1);
    if Xactive == 1, X = 1 + randn(T, 1); else, X = []; end
    D = zeros(T, 1); Y = zeros(T, 1);
    idx_r1 = floor(T/4);
    idx_r2 = idx_r1 + floor((1 - pi0) * T);
    for t = 1:T
        R_err = mvnrnd([0 0], Sigma_eu, 1);
        if t <= idx_r1 || t > idx_r2
            D(t) = theta_strong * Z(t) + R_err(1);
        else
            D(t) = theta_weak * Z(t) + R_err(1);
        end
        if Xactive == 1, D(t) = D(t) + gamma * X(t); end
        Y(t) = beta * D(t) + R_err(2);
        if Xactive == 1, Y(t) = Y(t) + gamma * X(t); end
    end
end


%% ========================================================================
%  EXPORT MODULE (EXACT PIVOTING TO PAPER FORMAT, EXTENDED WITH MM COLUMNS)
% =========================================================================

function export_paper_exhibits_sec6_MM_to_txt(Results, n_sims, folder_name, run_target, ts_file, ts_disp, rej_names_all, raw_names_all)
    filename = fullfile(folder_name, sprintf('Paper_Exhibits_SEC6_%dsims_MM_%s.txt', n_sims, ts_file));
    fid = fopen(filename, 'w');
    if fid == -1, return; end

    % Display labels: legacy MATLAB field name -> paper-style hyphenated label (see contract #9)
    display_labels = struct( ...
        'AR_AMS','AR_T', 'LM_AMS','LM_T', 'LR_AMS','CLR_T', ...
        'AR_ST','AR_T(\hat{S}_T)', 'LM_ST','LM_T(\hat{S}_T)', 'LR_ST','CLR_T(\hat{S}_T)', ...
        'GMM_S_MM','GMM-S', 'split_S_MM_Feasible','split-S (GMM)', 'qLL_S_MM','qLL-S', ...
        'exp_S_MM','exp-S', 'ave_S_MM','ave-S', ...
        'split_S_IV','split-S (IV)', 'split_KLM_IV','split-KLM', 'split_JKLM_IV','split-JKLM', 'split_CLR_IV','split-CLR', ...
        'M2_AMS','M2(AMS)', 'M2_ST','M2(S_hat)', 'rk_stat_IV','rk-stat' ...
    );

    fprintf(fid, '=======================================================================================================\n');
    fprintf(fid, '                 REPLICATION OF SECTION 6 PAPER EXHIBITS -- MM-EXTENDED (AMS/ST + Magnusson-Mavroeidis)\n');
    fprintf(fid, '=======================================================================================================\n');
    fprintf(fid, 'MC Iterations      : %d\n', n_sims);
    fprintf(fid, 'Date of Execution  : %s\n', ts_disp);
    fprintf(fid, 'Target Executed    : %s\n', upper(run_target));
    fprintf(fid, 'Rejection stats (15): %s\n', strjoin(rej_names_all, ', '));
    fprintf(fid, 'Raw diagnostics (3) : %s  [no rejection rule -- see contract #7 in runSection6_simulMM.m]\n', strjoin(raw_names_all, ', '));
    fprintf(fid, '=======================================================================================================\n\n');

    run_t9  = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab9');
    run_t10 = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab10');
    run_f7  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig7');
    run_f8  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig8');

    rej_labels = cellfun(@(f) display_labels.(f), rej_names_all, 'UniformOutput', false);
    raw_labels = cellfun(@(f) display_labels.(f), raw_names_all, 'UniformOutput', false);

    % --- TABLE 9 PIVOTING ---
    if run_t9
        fprintf(fid, 'Table 9 (MM-extended): Finite-Sample Null Rejection Frequencies of Tests\n');
        fprintf(fid, 'Model M1 and M2. The null hypothesis is H0: \\beta=0\n\n');
        pi0_grid = [0.4, 0.6, 0.8]; d_grid = [4, 8, 12, 16];
        for pi0 = pi0_grid
            fprintf(fid, '=======================================================================================================\n');
            fprintf(fid, 'T=400, \\rho=0.25, \\pi_0=%.1f\n', pi0);
            fprintf(fid, '-------------------------------------------------------------------------------------------------------\n');
            fprintf(fid, '%-18s | Model M1 (X0), d= %s | Model M2 (X1), d= %s\n', 'Statistic', strjoin(arrayfun(@(x) sprintf('%d',x), d_grid, 'UniformOutput', false), ' '), strjoin(arrayfun(@(x) sprintf('%d',x), d_grid, 'UniformOutput', false), ' '));
            fprintf(fid, '-------------------------------------------------------------------------------------------------------\n');
            for k = 1:length(rej_names_all)
                fprintf(fid, '%-18s', rej_labels{k});
                for d = d_grid
                    lbl = sprintf('Tab9_d%d_pi0%02d_X0', d, pi0*100);
                    if isfield(Results, lbl), val = Results.(lbl).size(k); else, val = NaN; end
                    fprintf(fid, ' %6.4f', val);
                end
                fprintf(fid, '  |');
                for d = d_grid
                    lbl = sprintf('Tab9_d%d_pi0%02d_X1', d, pi0*100);
                    if isfield(Results, lbl), val = Results.(lbl).size(k); else, val = NaN; end
                    fprintf(fid, ' %6.4f', val);
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
            fprintf(fid, 'Raw diagnostics (mean, no rejection rule):\n');
            for k = 1:length(raw_names_all)
                fprintf(fid, '%-18s', raw_labels{k});
                for d = d_grid
                    lbl = sprintf('Tab9_d%d_pi0%02d_X0', d, pi0*100);
                    if isfield(Results, lbl), val = Results.(lbl).raw(k); else, val = NaN; end
                    fprintf(fid, ' %6.4f', val);
                end
                fprintf(fid, '  |');
                for d = d_grid
                    lbl = sprintf('Tab9_d%d_pi0%02d_X1', d, pi0*100);
                    if isfield(Results, lbl), val = Results.(lbl).raw(k); else, val = NaN; end
                    fprintf(fid, ' %6.4f', val);
                end
                fprintf(fid, '\n');
            end
        end
        fprintf(fid, '=======================================================================================================\n\n');
    end

    % --- TABLE 10 PIVOTING ---
    if run_t10
        fprintf(fid, 'Table 10 (MM-extended): Finite-Sample Null Rejection Frequencies of Tests\n\n');
        X_grid = [0, 1]; rho_grid = [0.25, 0.50, 0.75]; d_grid = [4, 10, 16];
        for X = X_grid
            fprintf(fid, '=========================================================================================================================\n');
            if X == 0, fprintf(fid, '\\gamma_1=\\gamma_2=0 and \\rho_e=\\rho_u=0 (Model M1)\n'); else, fprintf(fid, '\\gamma_1=\\gamma_2=1 and \\rho_e=\\rho_u=0 (Model M2)\n'); end
            fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(fid, 'T=200, \\pi_0=0.6    \\rho=0.25 | \\rho=0.50 | \\rho=0.75   (d = 4 10 16 within each block)\n');
            fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------\n');
            for k = 1:length(rej_names_all)
                fprintf(fid, '%-18s', rej_labels{k});
                for rho = rho_grid
                    for d = d_grid
                        lbl = sprintf('Tab10_X%d_rho%02d_d%d', X, rho*100, d);
                        if isfield(Results, lbl), val = Results.(lbl).size(k); else, val = NaN; end
                        fprintf(fid, ' %6.4f', val);
                    end
                    if rho ~= 0.75, fprintf(fid, '  |'); end
                end
                fprintf(fid, '\n');
            end
            fprintf(fid, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n');
            fprintf(fid, 'Raw diagnostics (mean, no rejection rule):\n');
            for k = 1:length(raw_names_all)
                fprintf(fid, '%-18s', raw_labels{k});
                for rho = rho_grid
                    for d = d_grid
                        lbl = sprintf('Tab10_X%d_rho%02d_d%d', X, rho*100, d);
                        if isfield(Results, lbl), val = Results.(lbl).raw(k); else, val = NaN; end
                        fprintf(fid, ' %6.4f', val);
                    end
                    if rho ~= 0.75, fprintf(fid, '  |'); end
                end
                fprintf(fid, '\n');
            end
        end
        fprintf(fid, '=========================================================================================================================\n\n');
    end

    % --- FIGURE 7 PIVOTING ---
    fields_all = fieldnames(Results);
    if run_f7
        fig7_fields = fields_all(startsWith(fields_all, 'Fig7_'));
        if ~isempty(fig7_fields)
            fprintf(fid, '============================================================================================================\n');
            fprintf(fid, ' FIGURE 7 (MM-extended): Size-Adjusted Power Curve (Standard Power)\n');
            fprintf(fid, '============================================================================================================\n');
            for i = 1:length(fig7_fields)
                res = Results.(fig7_fields{i});
                fprintf(fid, '>>> Scenario: %s\n', fig7_fields{i});
                fprintf(fid, '      %-6s | %s\n', 'Beta', strjoin(rej_labels, ' | '));
                fprintf(fid, '      ----------------------------------------------------------------------------------------------\n');
                for b_idx = 1:length(res.beta_grid)
                    beta_lbl = sprintf('%4.2f', res.beta_grid(b_idx));
                    if res.beta_grid(b_idx) == 0.0, beta_lbl = 'SIZE'; end
                    row_str = strjoin(arrayfun(@(x) sprintf('%6.4f',x), res.power(b_idx,:), 'UniformOutput', false), ' | ');
                    fprintf(fid, '      %-6s | %s\n', beta_lbl, row_str);
                end
                fprintf(fid, '      -- Raw diagnostics (mean trajectory, no rejection rule: %s) --\n', strjoin(raw_labels, ', '));
                for b_idx = 1:length(res.beta_grid)
                    beta_lbl = sprintf('%4.2f', res.beta_grid(b_idx));
                    if res.beta_grid(b_idx) == 0.0, beta_lbl = 'SIZE'; end
                    row_str = strjoin(arrayfun(@(x) sprintf('%6.4f',x), res.raw(b_idx,:), 'UniformOutput', false), ' | ');
                    fprintf(fid, '      %-6s | %s\n', beta_lbl, row_str);
                end
                fprintf(fid, '\n');
            end
        end
    end

    % --- FIGURE 8 PIVOTING ---
    if run_f8
        fig8_fields = fields_all(startsWith(fields_all, 'Fig8_'));
        if ~isempty(fig8_fields)
            fprintf(fid, '============================================================================================================\n');
            fprintf(fid, ' FIGURE 8 (MM-extended): Size-Adjusted Power Curve (Perverse Identification)\n');
            fprintf(fid, '============================================================================================================\n');
            for i = 1:length(fig8_fields)
                res = Results.(fig8_fields{i});
                fprintf(fid, '>>> Scenario: %s\n', fig8_fields{i});
                fprintf(fid, '      %-6s | %s\n', 'Beta', strjoin(rej_labels, ' | '));
                fprintf(fid, '      ----------------------------------------------------------------------------------------------\n');
                for b_idx = 1:length(res.beta_grid)
                    beta_lbl = sprintf('%4.2f', res.beta_grid(b_idx));
                    if res.beta_grid(b_idx) == 0.0, beta_lbl = 'SIZE'; end
                    row_str = strjoin(arrayfun(@(x) sprintf('%6.4f',x), res.power(b_idx,:), 'UniformOutput', false), ' | ');
                    fprintf(fid, '      %-6s | %s\n', beta_lbl, row_str);
                end
                fprintf(fid, '      -- Raw diagnostics (mean trajectory, no rejection rule: %s) --\n', strjoin(raw_labels, ', '));
                for b_idx = 1:length(res.beta_grid)
                    beta_lbl = sprintf('%4.2f', res.beta_grid(b_idx));
                    if res.beta_grid(b_idx) == 0.0, beta_lbl = 'SIZE'; end
                    row_str = strjoin(arrayfun(@(x) sprintf('%6.4f',x), res.raw(b_idx,:), 'UniformOutput', false), ' | ');
                    fprintf(fid, '      %-6s | %s\n', beta_lbl, row_str);
                end
                fprintf(fid, '\n');
            end
        end
    end

    fclose(fid);
end
