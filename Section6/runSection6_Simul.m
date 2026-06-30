function Paper_Results = runSection6_Simul(n_sims_input, run_target)
% RUNSECTION6_SIMUL Executes the Monte Carlo suite for Section 6 exhibits.
% FOCUS: AMS and ST tests (pi-LATE methods). MM tests excluded.
%
% SYNTAX:
%   res = runSection6_Simul()                   % Default: All (10000 sims)
%   res = runSection6_Simul(10000, 'tab9')      % Runs EXCLUSIVELY Table 9
%   res = runSection6_Simul(10000, 'tab10')     % Runs EXCLUSIVELY Table 10
%   res = runSection6_Simul(10000, 'fig8')      % Runs EXCLUSIVELY Figure 8
%
% INPUTS:
%   n_sims_input : Number of Monte Carlo replications (Default: 10000)
%   run_target   : Target exhibit -> 'all', 'tab9', 'tab10', 'fig7', 'fig8'
% =========================================================================

    clc; 
    restoredefaultpath;
    addpath('..\m_break'); % Requires m_break in the MATLAB path
    warning('off','all');

    % --- OPTIONAL INPUTS AND DEFAULT HANDLING ---
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

    % --- DESCRIPTIVE STRINGS FOR LOGS AND DIRECTORY SETUP ---
    folder_name = sprintf('Paper_Exhibits_Results_%dsims_SEC6_%s', n_sims, upper(run_target));
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end

    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    timestamp_disp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    log_filename = fullfile(folder_name, sprintf('Detailed_Logs_%dsims_SEC6_%s_%s.txt', n_sims, upper(run_target), timestamp_file));
    fid_log = fopen(log_filename, 'w');
    if fid_log == -1, error('Unable to create the log file.'); end

    fprintf_both(fid_log, '===================================================================\n');
    fprintf_both(fid_log, ' STARTING SIMULATION SUITE: SECTION 6 (DYNAMIC LATE)\n');
    fprintf_both(fid_log, ' Execution Target: %s\n', upper(run_target));
    fprintf_both(fid_log, ' Date/Time       : %s\n', timestamp_disp);
    fprintf_both(fid_log, ' Iterations      : %d\n', n_sims);
    fprintf_both(fid_log, ' Focus           : AMS and ST Tests\n');
    fprintf_both(fid_log, ' Engine          : robust_LATE_inference (External Function)\n');
    fprintf_both(fid_log, '===================================================================\n\n');

    Paper_Results = struct();
    test_names = {'AR_AMS', 'LM_AMS', 'LR_AMS', 'AR_ST', 'LM_ST', 'LR_ST'};
    
    % Core logic parameters (legacy math framework structurally embedded)
    params_core = struct('q', 1, 'max_breaks', 5, 'eps1', 0.05, 'vector_piL', (0.6:0.1:1), 'step_pi', 0.05, 'Hac_option', true);

    % =========================================================================
    % ROUTINE: TABLE 9 (Empirical Size, i.i.d., T=400)
    % =========================================================================
    if run_tab9
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLE 9 (Empirical Size, i.i.d., T=400)\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        T_val = 400; beta_H0 = 0.0; rho_eu = 0.25;
        d_grid = [4, 8, 12, 16]; pi0_grid = [0.4, 0.6, 0.8]; Xactive_grid = [0, 1];
        
        for pi0_val = pi0_grid
            for Xactive = Xactive_grid
                for d_val = d_grid
                    run_label = sprintf('Tab9_d%d_pi0%02d_X%d', d_val, pi0_val*100, Xactive);
                    fprintf_both(fid_log, '\n--- Starting Iteration: Tab 9 [d=%d, pi0=%.1f, Xactive=%d] ---\n', d_val, pi0_val, Xactive);
                    
                    theta_H0 = d_val / sqrt(T_val);
                    rej_matrix = zeros(n_sims, 6);
                    
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = d_val*1e8 + (pi0_val*100)*1e6 + Xactive*1e5 + i_sims; 
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H0, 0, beta_H0, 1, rho_eu, Xactive, rng_seed);
                        [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        rej_matrix(i_sims, :) = [ ...
                            (stats_AMS.AR * params_core.q) > cv_LR.cv_AR, stats_AMS.LM > cv_LR.cv_LM, stats_AMS.LR > cv_LR.cv_AMS, ...
                            (stats_ST.AR * params_core.q)  > cv_LR.cv_AR, stats_ST.LM  > cv_LR.cv_LM, stats_ST.LR  > cv_LR.cv_ST ...
                        ];
                    end
                    emp_size = mean(rej_matrix, 1, 'omitnan');
                    Paper_Results.(run_label) = struct('params', [d_val, pi0_val, Xactive], 'size', emp_size);
                    
                    print_partial_results(fid_log, emp_size, test_names);
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
        
        for Xactive = Xactive_grid
            for rho_val = rho_grid
                for d_val = d_grid
                    run_label = sprintf('Tab10_X%d_rho%02d_d%d', Xactive, rho_val*100, d_val);
                    fprintf_both(fid_log, '\n--- Starting Iteration: Tab 10 [Xactive=%d, rho=%.2f, d=%d] ---\n', Xactive, rho_val, d_val);
                    
                    theta_H0 = d_val / sqrt(T_val);
                    rej_matrix = zeros(n_sims, 6);
                    
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = 1*1e6 + i_sims*100 + 16e7;
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H0, 0, beta_H0, 1, rho_val, Xactive, rng_seed);
                        [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        rej_matrix(i_sims, :) = [ ...
                            (stats_AMS.AR * params_core.q) > cv_LR.cv_AR, stats_AMS.LM > cv_LR.cv_LM, stats_AMS.LR > cv_LR.cv_AMS, ...
                            (stats_ST.AR * params_core.q)  > cv_LR.cv_AR, stats_ST.LM  > cv_LR.cv_LM, stats_ST.LR  > cv_LR.cv_ST ...
                        ];
                    end
                    emp_size = mean(rej_matrix, 1, 'omitnan');
                    Paper_Results.(run_label) = struct('params', [Xactive, rho_val, d_val], 'size', emp_size);
                    
                    print_partial_results(fid_log, emp_size, test_names);
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
        
        for rho_val = rho_grid
            for d_val = d_grid
                run_label = sprintf('Fig7_rho%02d_d%d', rho_val*100, d_val);
                fprintf_both(fid_log, '\n--- Starting Power Curve: Fig 7 [rho=%.2f, d=%d] ---\n', rho_val, d_val);
                
                theta_H1_strong = d_val / sqrt(T_val);
                stats_H0_matrix = zeros(n_sims, 6);
                
                fprintf_both(fid_log, '    > Computing empirical CVs under H0...\n');
                parfor i_sims = 1:n_sims
                    warning('off','all');
                    rng_seed = 1*1e6 + i_sims*100 + 16e7; 
                    
                    [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, 0, beta_H0, 1, rho_val, Xactive, rng_seed);
                    [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                    
                    stats_H0_matrix(i_sims, :) = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                end
                emp_CVs = prctile(stats_H0_matrix, 95, 1);
                
                power_matrix = zeros(length(beta_grid), 6);
                fprintf_both(fid_log, '    > Computing Power over beta grid...\n');
                for b_idx = 1:length(beta_grid)
                    beta_alt = beta_grid(b_idx);
                    rej_H1 = zeros(n_sims, 6);
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = b_idx*1e7 + i_sims*100 + 25e7; 
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, 0, beta_alt, 1, rho_val, Xactive, rng_seed);
                        [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                        rej_H1(i_sims, :) = vals > emp_CVs;
                    end
                    power_matrix(b_idx, :) = mean(rej_H1, 1, 'omitnan');
                end
                Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix);
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
        
        for p_idx = 1:length(perverse_types)
            p_type = perverse_types{p_idx};
            for pi0_val = pi0_grid
                for d_val = d_grid
                    
                    theta_H1_strong = d_val / sqrt(T_val);
                    if strcmp(p_type, 'local'), theta_H1_weak = -0.5 * d_val / sqrt(T_val);
                    else, theta_H1_weak = -0.5; end
                    
                    run_label = sprintf('Fig8_%s_pi0%02d_d%d', p_type, pi0_val*100, d_val);
                    fprintf_both(fid_log, '\n--- Starting Power Curve: Fig 8 [Type=%s, pi0=%.1f, d=%d] ---\n', p_type, pi0_val, d_val);
                    
                    stats_H0_matrix = zeros(n_sims, 6);
                    fprintf_both(fid_log, '    > Computing empirical CVs under H0...\n');
                    parfor i_sims = 1:n_sims
                        warning('off','all');
                        rng_seed = 1*1e6 + i_sims*100 + 16e7;
                        
                        [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, theta_H1_weak, beta_H0, 1, rho_val, Xactive, rng_seed);
                        [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                        
                        stats_H0_matrix(i_sims, :) = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                    end
                    emp_CVs = prctile(stats_H0_matrix, 95, 1);
                    
                    power_matrix = zeros(length(beta_grid), 6);
                    fprintf_both(fid_log, '    > Computing Power over beta grid...\n');
                    for b_idx = 1:length(beta_grid)
                        beta_alt = beta_grid(b_idx);
                        rej_H1 = zeros(n_sims, 6);
                        parfor i_sims = 1:n_sims
                            warning('off','all');
                            rng_seed = b_idx*1e7 + i_sims*100 + 25e7;
                            
                            [Y_s, D_s, Z_raw, X_raw] = generate_DGP(T_val, pi0_val, theta_H1_strong, theta_H1_weak, beta_alt, 1, rho_val, Xactive, rng_seed);
                            [stats_AMS, stats_ST, ~] = robust_LATE_inference(Y_s, D_s, Z_raw, X_raw, beta_H0, params_core, rng_seed);
                            
                            vals = [stats_AMS.AR, stats_AMS.LM, stats_AMS.LR, stats_ST.AR, stats_ST.LM, stats_ST.LR];
                            rej_H1(i_sims, :) = vals > emp_CVs;
                        end
                        power_matrix(b_idx, :) = mean(rej_H1, 1, 'omitnan');
                    end
                    Paper_Results.(run_label) = struct('beta_grid', beta_grid, 'power', power_matrix);
                    fprintf_both(fid_log, '    > Power computed successfully.\n');
                end
            end
        end
    end

    % --- FINAL SHUTDOWN AND EXPORT ---
    reset_parallel_pool(fid_log); 
    fclose(fid_log);
    
    out_filename_mat = fullfile(folder_name, sprintf('Paper_Exhibits_Suite_SEC6_%dsims_%s_%s.mat', n_sims, upper(run_target), timestamp_file));
    save(out_filename_mat, 'Paper_Results');
    
    export_paper_exhibits_sec6_to_txt(Paper_Results, n_sims, folder_name, run_target, timestamp_file, timestamp_disp);
    
    fprintf('\nExecution Completed Successfully (Target: %s). Outputs saved in: %s\\\n', upper(run_target), folder_name);
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