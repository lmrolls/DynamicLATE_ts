function Paper_Results = runSection6_Simul(n_sims_input, run_target, include_competitors)
% RUNSECTION6_SIMUL Executes the Monte Carlo suite for Section 6 exhibits.
% FOCUS: AMS and ST tests (pi-LATE methods), optionally together with the
% competitor tests of Magnusson and Mavroeidis (2014) ('MM', from competitors_MM\).
%
% SYNTAX:
%   res = runSection6_Simul()                   % Default: All (10000 sims)
%   res = runSection6_Simul(10000, 'tab9')      % Runs EXCLUSIVELY Table 9
%   res = runSection6_Simul(10000, 'tab10')     % Runs EXCLUSIVELY Table 10
%   res = runSection6_Simul(10000, 'fig8')      % Runs EXCLUSIVELY Figure 8
%   res = runSection6_Simul(10000, 'all', false)% pi-LATE tests only (no competitors)
%
% INPUTS:
%   n_sims_input : Number of Monte Carlo replications (Default: 10000)
%   run_target   : Target exhibit -> 'all', 'tab9', 'tab10', 'fig7', 'fig8'
%   include_competitors : (Optional) true (default) also computes the 9 MM test
%                  statistics (GMM-S, split-S (GMM), qLL-S, exp-S, ave-S, split-S (IV),
%                  split-KLM, split-JKLM, split-CLR) and the raw diagnostics
%                  (M2_AMS, M2_ST, rk-stat) with the legacy MM code in competitors_MM\;
%                  false reproduces exactly the previous pi-LATE-only outputs (same
%                  folder and file names, same fields, same tables).
%                  With true the output folder and file names carry the token _MM and
%                  'size' / 'power' have 15 columns (6 pi-LATE + 9 MM) plus 'raw' (3).
%
% COMPETITORS (only when include_competitors = true):
%   * MM statistics, HAC options and critical values exactly as in the legacy
%     scripts (NW kernel, lag floor(T^(1/3)), trimming 0.05, q=1, p=1, p_zeta=0).
%   * Asymptotic/simulated CVs of qLL-S, exp-S, ave-S: computed once per run;
%     for Table 9, as in the legacy script simulFun_TAB_9.m, they are re-simulated
%     for every scenario after rng(d*1000 + pi0*100 + Xactive) (set
%     mm_cv_seed_per_scenario = false below to keep the once-per-run values).
%   * split-CLR CV simulated per replication with seed rng_seed+3 (legacy convention).
%   * Figures 7-8: size-adjusted power of all statistics from the empirical 95th
%     percentile under H0, as for the pi-LATE tests.
%   * The pi-LATE stages, DGP and seeds are unchanged: the 6 pi-LATE columns are
%     identical with and without the competitors.
% =========================================================================

    clc; 
    restoredefaultpath;
    addpath('..\m_break'); % Requires m_break in the MATLAB path
    warning('off','all');

    % --- OPTIONAL INPUTS AND DEFAULT HANDLING ---
    if nargin < 1 || isempty(n_sims_input), n_sims = 10000; else, n_sims = n_sims_input; end
    if nargin < 2 || isempty(run_target), run_target = 'all'; end
    if nargin < 3 || isempty(include_competitors), include_competitors = true; end
    if include_competitors
        addpath('competitors_MM');  % computeMMtests2, computeMM_IVstats2, getMagnussonAsymptoticCV, split_CLR_cv_seed, HAC helpers
        comp_tag = '_MM';
    else
        comp_tag = '';
    end
    mm_cv_seed_per_scenario = true;   % Table 9: legacy per-scenario seeding of the MM asymptotic CVs
    
    valid_targets = {'all', 'tab9', 'tab10', 'fig7', 'fig8'};
    if ~ismember(lower(run_target), valid_targets)
        error('Error: run_target must be ''all'', ''tab9'', ''tab10'', ''fig7'', or ''fig8''.');
    end
    
    run_tab9  = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab9');
    run_tab10 = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab10');
    run_fig7  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig7');
    run_fig8  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig8');

    % --- DESCRIPTIVE STRINGS FOR LOGS AND DIRECTORY SETUP ---
    folder_name = sprintf('Paper_Exhibits_Results_%dsims_SEC6_%s%s', n_sims, upper(run_target), comp_tag);
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end

    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    timestamp_disp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    log_filename = fullfile(folder_name, sprintf('Detailed_Logs_%dsims_SEC6_%s%s_%s.txt', n_sims, upper(run_target), comp_tag, timestamp_file));
    fid_log = fopen(log_filename, 'w');
    if fid_log == -1, error('Unable to create the log file.'); end

    fprintf_both(fid_log, '===================================================================\n');
    if include_competitors
        fprintf_both(fid_log, ' STARTING SIMULATION SUITE: SECTION 6 (DYNAMIC LATE) -- MM-EXTENDED\n');
    else
        fprintf_both(fid_log, ' STARTING SIMULATION SUITE: SECTION 6 (DYNAMIC LATE)\n');
    end
    fprintf_both(fid_log, ' Execution Target: %s\n', upper(run_target));
    fprintf_both(fid_log, ' Date/Time       : %s\n', timestamp_disp);
    fprintf_both(fid_log, ' Iterations      : %d\n', n_sims);
    if include_competitors
        fprintf_both(fid_log, ' Focus           : AMS, ST (pi-LATE) AND Magnusson-Mavroeidis (MM) Tests\n');
        fprintf_both(fid_log, ' Engine          : robust_LATE_inference (unmodified) + computeMMtests2 / computeMM_IVstats2 (competitors_MM)\n');
    else
        fprintf_both(fid_log, ' Focus           : AMS and ST Tests\n');
        fprintf_both(fid_log, ' Engine          : robust_LATE_inference (External Function)\n');
    end
    fprintf_both(fid_log, '===================================================================\n\n');

    Paper_Results = struct();
    test_names = {'AR_AMS', 'LM_AMS', 'LR_AMS', 'AR_ST', 'LM_ST', 'LR_ST'};
    
    % Core logic parameters (legacy math framework structurally embedded)
    params_core = struct('q', 1, 'max_breaks', 5, 'eps1', 0.05, 'vector_piL', (0.6:0.1:1), 'step_pi', 0.05, 'Hac_option', true);

    % --- Competitors (MM): statistic names, constants and critical values ---
    n_stats = 6;
    rej_names_all = test_names; raw_names_all = {};
    alpha_test_MM = NaN; q_instr_MM = NaN; p_param_dim_MM = NaN; trim_prop_MM = NaN; theoretical_CVs_MM = struct();
    if include_competitors
        n_stats = 15;
        % --- Canonical ordered statistic names, used only for reporting/printing ---
        rej_names_all = {'AR_AMS','LM_AMS','LR_AMS','AR_ST','LM_ST','LR_ST', ...
                          'GMM_S_MM','split_S_MM_Feasible','qLL_S_MM','exp_S_MM','ave_S_MM', ...
                          'split_S_IV','split_KLM_IV','split_JKLM_IV','split_CLR_IV'};   % 15 columns, fixed order
        raw_names_all = {'M2_AMS','M2_ST','rk_stat_IV'};                                 % 3 columns, fixed order

        % =========================================================================
        % PRECOMPUTE MM ASYMPTOTIC / SIMULATED CRITICAL VALUES ONCE (legacy scripts: once per run)
        % =========================================================================
        alpha_test_MM   = 0.05;
        q_instr_MM      = 1;
        p_param_dim_MM  = 1;
        p_zeta_MM       = 0;
        N_cv_sims_MM    = 10000;
        T_cv_approx_MM  = 10000;
        trim_prop_MM    = 0.05;    % legacy value (not the 15% of the MM paper): kept as in the legacy scripts

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
        % split_CLR_IV has NO entry here: it is data-dependent, computed per replication via split_CLR_cv_seed

        fprintf_both(fid_log, '    GMM-S=%.4f  split-S(GMM)=%.4f  qLL-S=%.4f  exp-S=%.4f  ave-S=%.4f\n', ...
            theoretical_CVs_MM.GMM_S_MM, theoretical_CVs_MM.split_S_MM_Feasible, theoretical_CVs_MM.qLL_S_MM, theoretical_CVs_MM.exp_S_MM, theoretical_CVs_MM.ave_S_MM);
        fprintf_both(fid_log, '    split-S(IV)=%.4f  split-KLM=%.4f  split-JKLM=%.4f\n\n', ...
            theoretical_CVs_MM.split_S_IV, theoretical_CVs_MM.split_KLM_IV, theoretical_CVs_MM.split_JKLM_IV);
    end

    % =========================================================================
    % ROUTINE: TABLE 9 (Empirical Size, i.i.d., T=400)
    % =========================================================================
    if run_tab9
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLE 9 (Empirical Size, i.i.d., T=400)\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        T_val = 400; beta_H0 = 0.0; rho_eu = 0.25;
        d_grid = [4, 8, 12, 16]; pi0_grid = [0.4, 0.6, 0.8]; Xactive_grid = [0, 1];
        params_MM_base = struct(); params_MM_IV_base = struct();
        if include_competitors
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
        end
        
        for pi0_val = pi0_grid
            for Xactive = Xactive_grid
                for d_val = d_grid
                    run_label = sprintf('Tab9_d%d_pi0%02d_X%d', d_val, pi0_val*100, Xactive);
                    fprintf_both(fid_log, '\n--- Starting Iteration: Tab 9 [d=%d, pi0=%.1f, Xactive=%d] ---\n', d_val, pi0_val, Xactive);
                    
                    theta_H0 = d_val / sqrt(T_val);
                    rej_matrix = zeros(n_sims, n_stats);
                    raw_matrix = zeros(n_sims, 3);
                    if include_competitors && mm_cv_seed_per_scenario
                        % Legacy simulFun_TAB_9.m: rng(initial_seed) with initial_seed = d*1000 + pi0*100 + Xactive,
                        % then the asymptotic CVs of qLL-S, exp-S, ave-S are simulated for this scenario (same order).
                        rng(d_val*1000 + (pi0_val*100) + Xactive);
                        theoretical_CVs_MM.qLL_S_MM = getMagnussonAsymptoticCV('qLL-S', q_instr_MM, alpha_test_MM, N_cv_sims_MM, T_cv_approx_MM, trim_prop_MM, p_zeta_MM);
                        theoretical_CVs_MM.exp_S_MM = getMagnussonAsymptoticCV('exp-S', q_instr_MM, alpha_test_MM, N_cv_sims_MM, T_cv_approx_MM, trim_prop_MM, p_zeta_MM);
                        theoretical_CVs_MM.ave_S_MM = getMagnussonAsymptoticCV('ave-S', q_instr_MM, alpha_test_MM, N_cv_sims_MM, T_cv_approx_MM, trim_prop_MM, p_zeta_MM);
                        fprintf_both(fid_log, '    > MM asymptotic CVs for this scenario (seed %d): qLL-S=%.4f  exp-S=%.4f  ave-S=%.4f\n', ...
                            d_val*1000 + (pi0_val*100) + Xactive, theoretical_CVs_MM.qLL_S_MM, theoretical_CVs_MM.exp_S_MM, theoretical_CVs_MM.ave_S_MM);
                    end
                    
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = d_val*1e8 + (pi0_val*100)*1e6 + Xactive*1e5 + i_sims; 
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H0, 0, beta_H0, 1, rho_eu, Xactive, rng_seed);
                        [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        rej_row_1_6 = [ ...
                            (stats_AMS.AR * params_core.q) > cv_LR.cv_AR, stats_AMS.LM > cv_LR.cv_LM, stats_AMS.LR > cv_LR.cv_AMS, ...
                            (stats_ST.AR * params_core.q)  > cv_LR.cv_AR, stats_ST.LM  > cv_LR.cv_LM, stats_ST.LR  > cv_LR.cv_ST ...
                        ];
                        rej_row = rej_row_1_6; raw_row = [NaN, NaN, NaN];
                        if include_competitors
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
                            rej_row = [rej_row_1_6, rej_row_7_11, rej_row_12_15];
                            raw_row = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                        end
                        rej_matrix(i_sims, :) = rej_row;
                        raw_matrix(i_sims, :) = raw_row;
                    end
                    emp_size = mean(rej_matrix, 1, 'omitnan');
                    if include_competitors
                        raw_means = mean(raw_matrix, 1, 'omitnan');
                        Paper_Results.(run_label) = struct('params', [d_val, pi0_val, Xactive], 'size', emp_size, 'raw', raw_means);
                        print_partial_results_MM(fid_log, emp_size, raw_means, rej_names_all, raw_names_all);
                    else
                        Paper_Results.(run_label) = struct('params', [d_val, pi0_val, Xactive], 'size', emp_size);
                        print_partial_results(fid_log, emp_size, test_names);
                    end
                end
            end
        end
    end

    % =========================================================================
    % ROUTINE: TABLE 10 (Empirical Size, i.i.d., T=200)
    % =========================================================================
    if run_tab10
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLE 10 (Empirical Size, i.i.d. errors, T=200)\n');
        fprintf_both(fid_log, '>>> FORCING COMMON RANDOM NUMBERS (LEGACY SEED MATCH)\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        T_val = 200; beta_H0 = 0.0; pi0_val = 0.6;
        rho_grid = [0.25, 0.50, 0.75]; d_grid = [4, 10, 16]; Xactive_grid = [0, 1];
        params_MM_base = struct(); params_MM_IV_base = struct();
        if include_competitors
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
        end
        
        for Xactive = Xactive_grid
            for rho_val = rho_grid
                for d_val = d_grid
                    run_label = sprintf('Tab10_X%d_rho%02d_d%d', Xactive, rho_val*100, d_val);
                    fprintf_both(fid_log, '\n--- Starting Iteration: Tab 10 [Xactive=%d, rho=%.2f, d=%d] ---\n', Xactive, rho_val, d_val);
                    
                    theta_H0 = d_val / sqrt(T_val);
                    rej_matrix = zeros(n_sims, n_stats);
                    raw_matrix = zeros(n_sims, 3);
                    
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = 1*1e6 + i_sims*100 + 16e7;
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H0, 0, beta_H0, 1, rho_val, Xactive, rng_seed);
                        [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        rej_row_1_6 = [ ...
                            (stats_AMS.AR * params_core.q) > cv_LR.cv_AR, stats_AMS.LM > cv_LR.cv_LM, stats_AMS.LR > cv_LR.cv_AMS, ...
                            (stats_ST.AR * params_core.q)  > cv_LR.cv_AR, stats_ST.LM  > cv_LR.cv_LM, stats_ST.LR  > cv_LR.cv_ST ...
                        ];
                        rej_row = rej_row_1_6; raw_row = [NaN, NaN, NaN];
                        if include_competitors
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
                            rej_row = [rej_row_1_6, rej_row_7_11, rej_row_12_15];
                            raw_row = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                        end
                        rej_matrix(i_sims, :) = rej_row;
                        raw_matrix(i_sims, :) = raw_row;
                    end
                    emp_size = mean(rej_matrix, 1, 'omitnan');
                    if include_competitors
                        raw_means = mean(raw_matrix, 1, 'omitnan');
                        Paper_Results.(run_label) = struct('params', [Xactive, rho_val, d_val], 'size', emp_size, 'raw', raw_means);
                        print_partial_results_MM(fid_log, emp_size, raw_means, rej_names_all, raw_names_all);
                    else
                        Paper_Results.(run_label) = struct('params', [Xactive, rho_val, d_val], 'size', emp_size);
                        print_partial_results(fid_log, emp_size, test_names);
                    end
                end
            end
        end
    end

    % =========================================================================
    % ROUTINE: FIGURE 7 (Size-Adjusted Power Curve i.i.d. T=200)
    % =========================================================================
    if run_fig7
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: FIGURE 7 (Size-Adjusted Power, T=200)\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        T_val = 200; beta_H0 = 0.0; pi0_val = 0.6; Xactive = 0;
        rho_grid = [0.25, 0.75]; d_grid = [4, 10, 16]; beta_grid = 0:0.1:1;
        params_MM_base = struct(); params_MM_IV_base = struct();
        if include_competitors
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
        end
        
        for rho_val = rho_grid
            for d_val = d_grid
                run_label = sprintf('Fig7_rho%02d_d%d', rho_val*100, d_val);
                fprintf_both(fid_log, '\n--- Starting Power Curve: Fig 7 [rho=%.2f, d=%d] ---\n', rho_val, d_val);
                
                theta_H1_strong = d_val / sqrt(T_val);
                stats_H0_matrix = zeros(n_sims, n_stats);
                
                fprintf_both(fid_log, '    > Computing empirical CVs under H0...\n');
                parfor i_sims = 1:n_sims
                    warning('off','all');
                    rng_seed = 1*1e6 + i_sims*100 + 16e7; 
                    
                    [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, 0, beta_H0, 1, rho_val, Xactive, rng_seed);
                    [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                    
                    vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                    raw_row = [NaN, NaN, NaN];
                    if include_competitors
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
                        vals = [vals, ...
                            res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                            res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                        raw_row = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                    end
                    stats_H0_matrix(i_sims, :) = vals;
                end
                emp_CVs = prctile(stats_H0_matrix, 95, 1);
                
                power_matrix = zeros(length(beta_grid), n_stats);
                raw_matrix_fig = zeros(length(beta_grid), 3);
                fprintf_both(fid_log, '    > Computing Power over beta grid...\n');
                for b_idx = 1:length(beta_grid)
                    beta_alt = beta_grid(b_idx);
                    rej_H1 = zeros(n_sims, n_stats);
                    raw_H1 = zeros(n_sims, 3);
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = b_idx*1e7 + i_sims*100 + 25e7; 
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, 0, beta_alt, 1, rho_val, Xactive, rng_seed);
                        [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                        raw_row = [NaN, NaN, NaN];
                        if include_competitors
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
                            vals = [vals, ...
                                res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                                res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                            raw_row = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                        end
                        rej_H1(i_sims, :) = vals > emp_CVs;
                        raw_H1(i_sims, :) = raw_row;
                    end
                    power_matrix(b_idx, :) = mean(rej_H1, 1, 'omitnan');
                    raw_matrix_fig(b_idx, :) = mean(raw_H1, 1, 'omitnan');
                end
                if include_competitors
                    Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix, 'raw', raw_matrix_fig);
                else
                    Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix);
                end
                fprintf_both(fid_log, '    > Power computed successfully.\n');
            end
        end
    end

    % =========================================================================
    % ROUTINE: FIGURE 8 (Perverse Identification, i.i.d. T=200)
    % =========================================================================
    if run_fig8
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: FIGURE 8 (Perverse Identification, T=200)\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        T_val = 200; beta_H0 = 0.0; rho_val = 0.25; Xactive = 0; 
        pi0_grid = [0.6, 0.8]; d_grid = [24, 32]; beta_grid = 0:0.1:1;
        perverse_types = {'local', 'fixed'}; 
        params_MM_base = struct(); params_MM_IV_base = struct();
        if include_competitors
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
        end
        
        for p_idx = 1:length(perverse_types)
            p_type = perverse_types{p_idx};
            for pi0_val = pi0_grid
                for d_val = d_grid
                    
                    theta_H1_strong = d_val / sqrt(T_val);
                    if strcmp(p_type, 'local'), theta_H1_weak = -0.5 * d_val / sqrt(T_val);
                    else, theta_H1_weak = -0.5; end
                    
                    run_label = sprintf('Fig8_%s_pi0%02d_d%d', p_type, pi0_val*100, d_val);
                    fprintf_both(fid_log, '\n--- Starting Power Curve: Fig 8 [Type=%s, pi0=%.1f, d=%d] ---\n', p_type, pi0_val, d_val);
                    
                    stats_H0_matrix = zeros(n_sims, n_stats);
                    fprintf_both(fid_log, '    > Computing empirical CVs under H0...\n');
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = 1*1e6 + i_sims*100 + 16e7;
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, theta_H1_weak, beta_H0, 1, rho_val, Xactive, rng_seed);
                        [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                        raw_row = [NaN, NaN, NaN];
                        if include_competitors
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
                            vals = [vals, ...
                                res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                                res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                            raw_row = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                        end
                        stats_H0_matrix(i_sims, :) = vals;
                    end
                    emp_CVs = prctile(stats_H0_matrix, 95, 1);
                    
                    power_matrix = zeros(length(beta_grid), n_stats);
                    raw_matrix_fig = zeros(length(beta_grid), 3);
                    fprintf_both(fid_log, '    > Computing Power over beta grid...\n');
                    for b_idx = 1:length(beta_grid)
                        beta_alt = beta_grid(b_idx);
                        rej_H1 = zeros(n_sims, n_stats);
                        raw_H1 = zeros(n_sims, 3);
                        parfor i_sims = 1:n_sims
                            warning('off','all');
                            rng_seed = b_idx*1e7 + i_sims*100 + 25e7;
                            
                            [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, theta_H1_weak, beta_alt, 1, rho_val, Xactive, rng_seed);
                            [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                            
                            vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                            raw_row = [NaN, NaN, NaN];
                            if include_competitors
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
                                vals = [vals, ...
                                    res_MM.GMM_S_MM, res_MM.split_S_MM_Feasible, res_MM.qLL_S_MM, res_MM.exp_S_MM, res_MM.ave_S_MM, ...
                                    res_MM_IV.split_S_IV, res_MM_IV.split_KLM_IV, res_MM_IV.split_JKLM_IV, res_MM_IV.split_CLR_IV];
                                raw_row = [stats_AMS.M2, stats_ST.M2, res_MM_IV.rk_stat_IV];
                            end
                            rej_H1(i_sims, :) = vals > emp_CVs;
                            raw_H1(i_sims, :) = raw_row;
                        end
                        power_matrix(b_idx, :) = mean(rej_H1, 1, 'omitnan');
                        raw_matrix_fig(b_idx, :) = mean(raw_H1, 1, 'omitnan');
                    end
                    if include_competitors
                        Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix, 'raw', raw_matrix_fig);
                    else
                        Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix);
                    end
                    fprintf_both(fid_log, '    > Power computed successfully.\n');
                end
            end
        end
    end

    % --- FINAL SHUTDOWN AND EXPORT ---
    reset_parallel_pool(fid_log); 
    fclose(fid_log);
    
    out_filename_mat = fullfile(folder_name, sprintf('Paper_Exhibits_Suite_SEC6_%dsims_%s%s_%s.mat', n_sims, upper(run_target), comp_tag, timestamp_file));
    save(out_filename_mat, 'Paper_Results');
    
    if include_competitors
        export_paper_exhibits_sec6_MM_to_txt(Paper_Results, n_sims, folder_name, run_target, timestamp_file, timestamp_disp, rej_names_all, raw_names_all);
        fprintf('\nExecution Completed Successfully (Target: %s, MM-extended). Outputs saved in: %s\\\n', upper(run_target), folder_name);
    else
        export_paper_exhibits_sec6_to_txt(Paper_Results, n_sims, folder_name, run_target, timestamp_file, timestamp_disp);
        fprintf('\nExecution Completed Successfully (Target: %s). Outputs saved in: %s\\\n', upper(run_target), folder_name);
    end
end


%% ========================================================================
%  LOCAL UTILITIES (Data Generation)
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

function print_partial_results(fid, emp_size, test_names)
    fprintf_both(fid, '    > Partial Rejection Rates:\n');
    fprintf_both(fid, '      | %8s | %8s | %8s || %8s | %8s | %8s |\n', test_names{:});
    fprintf_both(fid, '      | %8.4f | %8.4f | %8.4f || %8.4f | %8.4f | %8.4f |\n\n', emp_size);
end

function print_partial_results_MM(fid, emp_size, raw_means, rej_names, raw_names)
    fprintf_both(fid, '    > Partial Rejection Rates (AMS/ST + MM, 15 statistics):\n');
    fprintf_both(fid, '      %s\n', strjoin(rej_names, ' | '));
    fprintf_both(fid, '      %s\n', strjoin(arrayfun(@(x) sprintf('%.4f',x), emp_size, 'UniformOutput', false), ' | '));
    fprintf_both(fid, '    > Raw diagnostic means (no rejection rule, %s):\n', strjoin(raw_names, ', '));
    fprintf_both(fid, '      %s\n\n', strjoin(arrayfun(@(x) sprintf('%.4f',x), raw_means, 'UniformOutput', false), ' | '));
end


function [Y, D, Z, X] = generate_DGP(T, pi0, theta_strong, theta_weak, beta, gamma, rho, Xactive, rng_seed)
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
%  EXPORT MODULE (EXACT PIVOTING TO PAPER FORMAT)
% =========================================================================

function export_paper_exhibits_sec6_to_txt(Results, n_sims, folder_name, run_target, ts_file, ts_disp)
    filename = fullfile(folder_name, sprintf('Paper_Exhibits_SEC6_%dsims_%s.txt', n_sims, ts_file));
    fid = fopen(filename, 'w');
    if fid == -1, return; end
    
    fprintf(fid, '=======================================================================================================\n');
    fprintf(fid, '                 REPLICATION OF SECTION 6 PAPER EXHIBITS\n');
    fprintf(fid, '=======================================================================================================\n');
    fprintf(fid, 'MC Iterations      : %d\n', n_sims);
    fprintf(fid, 'Date of Execution  : %s\n', ts_disp);
    fprintf(fid, 'Target Executed    : %s\n', upper(run_target));
    fprintf(fid, '=======================================================================================================\n\n');
    
    run_t9  = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab9');
    run_t10 = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab10');
    run_f7  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig7');
    run_f8  = strcmpi(run_target, 'all') || strcmpi(run_target, 'fig8');
    
    % Index mapping [AR_AMS, LM_AMS, LR_AMS, AR_ST, LM_ST, LR_ST]
    r_idx = [2, 3, 5, 6]; 
    r_names = {'LM_T', 'CLR_T', 'LM_T(\\hat{S}_T)', 'CLR_T(\\hat{S}_T)'};

    % --- TABLE 9 PIVOTING ---
    if run_t9
        fprintf(fid, 'Table 9: Finite-Sample Null Rejection Frequencies of Tests\n');
        fprintf(fid, 'Model M1 and M2. The null hypothesis is H0: \\beta=0\n\n');
        pi0_grid = [0.4, 0.6, 0.8]; d_grid = [4, 8, 12, 16];
        for pi0 = pi0_grid
            fprintf(fid, '=======================================================================================================\n');
            fprintf(fid, 'T=400, \\rho=0.25, \\pi_0=%.1f\n', pi0);
            fprintf(fid, '-------------------------------------------------------------------------------------------------------\n');
            fprintf(fid, '                    \\gamma_1=\\gamma_2=0 and \\rho_e=\\rho_u=0          | \\gamma_1=\\gamma_2=1 and \\rho_e=\\rho_u=0\n');
            fprintf(fid, 'Test                d=4     d=8     d=12    d=16    | d=4     d=8     d=12    d=16\n');
            fprintf(fid, '-------------------------------------------------------------------------------------------------------\n');
            for k = 1:length(r_idx)
                fprintf(fid, '%-19s', r_names{k});
                for d = d_grid
                    lbl = sprintf('Tab9_d%d_pi0%02d_X0', d, pi0*100);
                    if isfield(Results, lbl), val = Results.(lbl).size(r_idx(k)); else, val = NaN; end
                    fprintf(fid, ' %5.3f  ', val);
                end
                fprintf(fid, ' |');
                for d = d_grid
                    lbl = sprintf('Tab9_d%d_pi0%02d_X1', d, pi0*100);
                    if isfield(Results, lbl), val = Results.(lbl).size(r_idx(k)); else, val = NaN; end
                    fprintf(fid, ' %5.3f  ', val);
                end
                fprintf(fid, '\n');
            end
        end
        fprintf(fid, '=======================================================================================================\n\n');
    end

    % --- TABLE 10 PIVOTING ---
    if run_t10
        fprintf(fid, 'Table 10: Finite-Sample Null Rejection Frequencies of Tests\n\n');
        X_grid = [0, 1]; rho_grid = [0.25, 0.50, 0.75]; d_grid = [4, 10, 16];
        for X = X_grid
            fprintf(fid, '=========================================================================================================================\n');
            if X == 0, fprintf(fid, '\\gamma_1=\\gamma_2=0 and \\rho_e=\\rho_u=0\n'); else, fprintf(fid, '\\gamma_1=\\gamma_2=1 and \\rho_e=\\rho_u=0\n'); end
            fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(fid, 'T=200, \\pi_0=0.6    \\rho=0.25                       | \\rho=0.50                       | \\rho=0.75\n');
            fprintf(fid, 'Test                d=4     d=10    d=16    | d=4     d=10    d=16    | d=4     d=10    d=16\n');
            fprintf(fid, '-------------------------------------------------------------------------------------------------------------------------\n');
            for k = 1:length(r_idx)
                fprintf(fid, '%-19s', r_names{k});
                for rho = rho_grid
                    for d = d_grid
                        lbl = sprintf('Tab10_X%d_rho%02d_d%d', X, rho*100, d);
                        if isfield(Results, lbl), val = Results.(lbl).size(r_idx(k)); else, val = NaN; end
                        fprintf(fid, ' %5.3f  ', val);
                    end
                    if rho ~= 0.75, fprintf(fid, ' |'); end
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
            fprintf(fid, ' FIGURE 7: Size-Adjusted Power Curve (Standard Power)\n');
            fprintf(fid, '============================================================================================================\n');
            for i = 1:length(fig7_fields)
                res = Results.(fig7_fields{i});
                fprintf(fid, '>>> Scenario: %s\n', fig7_fields{i});
                fprintf(fid, '      %4s | %6s | %6s | %6s || %6s | %6s | %6s |\n', 'Beta', 'AR_AMS', 'LM_AMS', 'LR_AMS', 'AR_ST', 'LM_ST', 'LR_ST');
                fprintf(fid, '      -----------------------------------------------------------------------\n');
                for b_idx = 1:length(res.beta_grid)
                    beta_lbl = sprintf('%4.2f', res.beta_grid(b_idx));
                    if res.beta_grid(b_idx) == 0.0, beta_lbl = 'SIZE'; end
                    fprintf(fid, '      %4s | %6.4f | %6.4f | %6.4f || %6.4f | %6.4f | %6.4f |\n', ...
                        beta_lbl, res.power(b_idx,1), res.power(b_idx,2), res.power(b_idx,3), ...
                        res.power(b_idx,4), res.power(b_idx,5), res.power(b_idx,6));
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
            fprintf(fid, ' FIGURE 8: Size-Adjusted Power Curve (Perverse Identification)\n');
            fprintf(fid, '============================================================================================================\n');
            for i = 1:length(fig8_fields)
                res = Results.(fig8_fields{i});
                fprintf(fid, '>>> Scenario: %s\n', fig8_fields{i});
                fprintf(fid, '      %4s | %6s | %6s | %6s || %6s | %6s | %6s |\n', 'Beta', 'AR_AMS', 'LM_AMS', 'LR_AMS', 'AR_ST', 'LM_ST', 'LR_ST');
                fprintf(fid, '      -----------------------------------------------------------------------\n');
                for b_idx = 1:length(res.beta_grid)
                    beta_lbl = sprintf('%4.2f', res.beta_grid(b_idx));
                    if res.beta_grid(b_idx) == 0.0, beta_lbl = 'SIZE'; end
                    fprintf(fid, '      %4s | %6.4f | %6.4f | %6.4f || %6.4f | %6.4f | %6.4f |\n', ...
                        beta_lbl, res.power(b_idx,1), res.power(b_idx,2), res.power(b_idx,3), ...
                        res.power(b_idx,4), res.power(b_idx,5), res.power(b_idx,6));
                end
                fprintf(fid, '\n');
            end
        end
    end

    fclose(fid);
end


%% ========================================================================
%  EXPORT MODULE (MM-EXTENDED: 6 pi-LATE + 9 MM COLUMNS, RAW DIAGNOSTICS)
% =========================================================================

function export_paper_exhibits_sec6_MM_to_txt(Results, n_sims, folder_name, run_target, ts_file, ts_disp, rej_names_all, raw_names_all)
    filename = fullfile(folder_name, sprintf('Paper_Exhibits_SEC6_%dsims_MM_%s.txt', n_sims, ts_file));
    fid = fopen(filename, 'w');
    if fid == -1, return; end

    % Display labels: legacy MATLAB field name -> paper-style hyphenated label
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
    fprintf(fid, 'Raw diagnostics (3) : %s  [no rejection rule]\n', strjoin(raw_names_all, ', '));
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