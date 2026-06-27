function Paper_Results = runSection5_Simul(n_sims_input, run_target)
% RUNSECTION5_SIMUL Executes the Monte Carlo suite for Section 5 exhibits.
%
% SYNTAX:
%   res = runSection5_Simul()                   % Default: All tables (5000 sims)
%   res = runSection5_Simul(5000, 'tab67')      % Runs EXCLUSIVELY Tables 6 and 7
%   res = runSection5_Simul(5000, 'tab8')       % Runs EXCLUSIVELY Table 8
%
% INPUTS:
%   n_sims_input : Number of Monte Carlo replications (Default: 5000)
%   run_target   : Target exhibit -> 'all' (default), 'tab67', 'tab8'
%
% OUTPUTS:
%   Paper_Results: Struct containing all Bias and MSE metrics for the grid.
% =========================================================================

    clc;
    restoredefaultpath;
    warning('off','all');

    % --- DEFAULT INPUT HANDLING ---
    if nargin < 1 || isempty(n_sims_input), n_sims = 5000; else, n_sims = n_sims_input; end
    if nargin < 2 || isempty(run_target), run_target = 'all'; end
    
    valid_targets = {'all', 'tab67', 'tab8'};
    if ~ismember(lower(run_target), valid_targets)
        error('Error: run_target must be ''all'', ''tab67'', or ''tab8''.');
    end
    
    % Routing booleans
    run_tab67 = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab67');
    run_tab8  = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab8');

    % --- UNDATED DIRECTORY & DATED FILES SETUP ---
    folder_name = sprintf('Paper_Exhibits_Results_%dsims_SEC5_%s', n_sims, upper(run_target));
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end

    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    timestamp_disp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    log_filename = fullfile(folder_name, sprintf('Detailed_Logs_%dsims_SEC5_%s_%s.txt', n_sims, upper(run_target), timestamp_file));
    fid_log = fopen(log_filename, 'w');
    if fid_log == -1, error('Unable to create the log file.'); end

    fprintf_both(fid_log, '===================================================================\n');
    fprintf_both(fid_log, ' STARTING SIMULATION SUITE: SECTION 5 (ESTIMATION)\n');
    fprintf_both(fid_log, ' Execution Target: %s\n', upper(run_target));
    fprintf_both(fid_log, ' Date/Time       : %s\n', timestamp_disp);
    fprintf_both(fid_log, ' Iterations      : %d\n', n_sims);
    fprintf_both(fid_log, '===================================================================\n\n');

    Paper_Results = struct();
    T = 200;

    % =========================================================================
    % ROUTINE: TABLES 6 & 7
    % =========================================================================
    if run_tab67
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLES 6 & 7\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        pi0_vals = [0.6, 0.8];
        rho_vals = [0.25, 0.75];
        rho_eu_vals = [0.0, 0.5];
        
        for p = 1:length(pi0_vals)
            pi0 = pi0_vals(p);
            for r = 1:length(rho_vals)
                rho = rho_vals(r);
                for eu = 1:length(rho_eu_vals)
                    rho_eu = rho_eu_vals(eu);
                    
                    run_label = sprintf('Tab67_rho%02d_pi%02d_eu%02d', round(rho*100), round(pi0*100), round(rho_eu*100));
                    fprintf_both(fid_log, '\n--- Starting Iteration: %s [rho=%.2f, pi0=%.1f, rho_eu=%.1f] ---\n', run_label, rho, pi0, rho_eu);
                    
                    res = Run_Sec5_Scenario('tab67', rho, pi0, rho_eu, T, n_sims, fid_log);
                    Paper_Results.(run_label) = res;
                end
            end
        end
    end

    % =========================================================================
    % ROUTINE: TABLE 8
    % =========================================================================
    if run_tab8
        reset_parallel_pool(fid_log);
        fprintf_both(fid_log, '\n>>> INITIALIZING REPLICATION: TABLE 8\n');
        fprintf_both(fid_log, '===================================================================\n');
        
        dgp_types = {'tab8_upper', 'tab8_lower'};
        pi0 = 0.6;
        rho_vals = [0.25, 0.75];
        rho_eu = 0.0;
        
        for dt = 1:length(dgp_types)
            dgp_type = dgp_types{dt};
            if strcmp(dgp_type, 'tab8_upper')
                target_str = 'upper';
            else
                target_str = 'lower';
            end
            
            for r = 1:length(rho_vals)
                rho = rho_vals(r);
                run_label = sprintf('Tab8_%s_rho%02d', target_str, round(rho*100));
                fprintf_both(fid_log, '\n--- Starting Iteration: %s [DGP=%s, rho=%.2f] ---\n', run_label, dgp_type, rho);
                
                res = Run_Sec5_Scenario(dgp_type, rho, pi0, rho_eu, T, n_sims, fid_log);
                Paper_Results.(run_label) = res;
            end
        end
    end

    % --- SHUTDOWN AND SAVE ---
    reset_parallel_pool(fid_log);
    fclose(fid_log);
    
    % Save .MAT workspace
    out_filename_mat = fullfile(folder_name, sprintf('Paper_Exhibits_Suite_SEC5_%dsims_%s_%s.mat', n_sims, upper(run_target), timestamp_file));
    save(out_filename_mat, 'Paper_Results');
    
    % Export formatted LaTeX-style TXT table
    export_paper_exhibits_sec5_to_txt(Paper_Results, n_sims, folder_name, run_target, timestamp_file, timestamp_disp);
    
    fprintf('\nExecution Completed Successfully (Target: %s). Outputs saved in: %s\n', upper(run_target), folder_name);
end

%% ========================================================================
%  LOCAL FUNCTIONS (ORCHESTRATOR WRAPPER)
% =========================================================================

function reset_parallel_pool(fid)
    fprintf_both(fid, '\n[SYSTEM] Executing memory cleanup and resetting Parallel Pool...\n');
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
        fprintf_both(fid, '[SYSTEM] Previous Parallel Pool closed successfully.\n');
    else
        fprintf_both(fid, '[SYSTEM] No active Parallel Pool detected.\n');
    end
    total_cores = feature('numcores');
    num_workers = max(1, total_cores - 2); 
    parpool('local', num_workers);
    fprintf_both(fid, '[SYSTEM] Parallel Pool started with %d workers.\n', num_workers);
end

function fprintf_both(fid, fmt, varargin)
    fprintf(1, fmt, varargin{:});
    if fid ~= -1, fprintf(fid, fmt, varargin{:}); end
end

%% ========================================================================
%  CORE SIMULATION ENGINE (Encapsulating 3-parfor logic)
% =========================================================================

function res = Run_Sec5_Scenario(dgp_type, rho, pi0, rho_eu, T, numMCreps, fid)
    
    % Fixed Parameters (Verbatim from legacy script)
    gamma1 = 1; gamma2 = 1; X_active = 1; fixed_m = 2; burn_in = 100;
    beta_dgp = 1; beta_idx = 1;
    
    if strcmp(dgp_type, 'tab67')
        d_grid = [4, 16, 32]; d2_val = 0; 
    elseif strcmp(dgp_type, 'tab8_upper')
        d_grid = [16, 24, 32]; d2_val = 4;
    elseif strcmp(dgp_type, 'tab8_lower')
        d_grid = [16, 24, 32]; d2_val = 8;
    else
        error('Unrecognized DGP type.');
    end
    
    num_d = length(d_grid);
    Bias_OLS = zeros(1, num_d); MSE_OLS = zeros(1, num_d);
    Bias_FGLS = zeros(1, num_d); MSE_FGLS = zeros(1, num_d);
    MSE_S_OLS = zeros(1, num_d); MSE_S_FGLS = zeros(1, num_d);
    
    for d_idx = 1:num_d
        d_val = d_grid(d_idx);
        fprintf_both(fid, '  > Running d_val = %d ...\n', d_val);
        
        Y_all = cell(numMCreps, 1); D_all = cell(numMCreps, 1);
        Z_all = cell(numMCreps, 1); X_all = cell(numMCreps, 1);
        Index_S0_all = cell(numMCreps, 1);
        
        b_ols = zeros(numMCreps, 1);
        b_fgls = zeros(numMCreps, 1);
        mse_s_ols = zeros(numMCreps, 1);
        mse_s_fgls = zeros(numMCreps, 1);
        
        % ---------------------------------------------------------
        % STAGE 1: Data Generating Process
        % ---------------------------------------------------------
        parfor i_sim = 1:numMCreps
            if strcmp(dgp_type, 'tab67')
                [Y_all{i_sim}, D_all{i_sim}, Z_all{i_sim}, X_all{i_sim}, Index_S0_all{i_sim}] = ...
                    generate_data_tab67(T, pi0, beta_dgp, d_val, rho, rho_eu, rho_eu, burn_in, X_active, gamma1, gamma2, d_idx, beta_idx, i_sim);
            else
                [Y_all{i_sim}, D_all{i_sim}, Z_all{i_sim}, X_all{i_sim}, Index_S0_all{i_sim}] = ...
                    generate_data_tab8(T, pi0, beta_dgp, d_val, d2_val, rho, X_active, gamma1, gamma2, d_idx, beta_idx, i_sim);
            end
        end
        
        % ---------------------------------------------------------
        % STAGE 2: 2SLS pi-LATE 
        % ---------------------------------------------------------
        addpath('..\m_break');
        parfor i_sim = 1:numMCreps
            Y_sim = Y_all{i_sim}; D_sim = D_all{i_sim}; 
            Z_sim = Z_all{i_sim}; X_sim = X_all{i_sim}; 
            S0_sim = Index_S0_all{i_sim};
            
            % Greedy Optimal 2SLS
            [b_o_g, s_o_g, ~] = est_2sls_pilate(Y_sim, D_sim, Z_sim, X_sim, pi0, fixed_m);
            
            b_ols(i_sim) = b_o_g;
            mse_s_ols(i_sim) = sum((s_o_g - S0_sim).^2) / T;
        end
        rmpath('..\m_break');
        
        % ---------------------------------------------------------
        % STAGE 3: Multivariate FGLS pi-LATE
        % ---------------------------------------------------------
        addpath('..\Mbreak_Multivariate');
        parfor i_sim = 1:numMCreps
            Y_sim = Y_all{i_sim}; D_sim = D_all{i_sim}; 
            Z_sim = Z_all{i_sim}; X_sim = X_all{i_sim}; 
            S0_sim = Index_S0_all{i_sim};
            
            [b_f_g_d, s_f_g_d, ~] = est_fgls_pilate(Y_sim, D_sim, Z_sim, X_sim, pi0, fixed_m);
            
            b_fgls(i_sim) = b_f_g_d;
            mse_s_fgls(i_sim) = sum((s_f_g_d - S0_sim).^2) / T;
        end
        rmpath('..\Mbreak_Multivariate');
        
        % ---------------------------------------------------------
        % AGGREGATION & METRICS
        % ---------------------------------------------------------
        Bias_OLS(d_idx) = mean(b_ols - beta_dgp, 'omitnan');
        MSE_OLS(d_idx) = mean((b_ols - beta_dgp).^2, 'omitnan');
        
        Bias_FGLS(d_idx) = mean(b_fgls - beta_dgp, 'omitnan');
        MSE_FGLS(d_idx) = mean((b_fgls - beta_dgp).^2, 'omitnan');
        
        MSE_S_OLS(d_idx) = mean(mse_s_ols, 'omitnan');
        MSE_S_FGLS(d_idx) = mean(mse_s_fgls, 'omitnan');
        
        % ---------------------------------------------------------
        % VERBOSE LOGGING (Console & TXT)
        % ---------------------------------------------------------
        fprintf_both(fid, '    >> Completed d = %d. Partial Results (x 10^2):\n', d_val);
        fprintf_both(fid, '       * MSE(S_hat)  | OLS: %8.2f | FGLS: %8.2f\n', MSE_S_OLS(d_idx)*100, MSE_S_FGLS(d_idx)*100);
        fprintf_both(fid, '       * Bias(beta)  | OLS: %8.2f | FGLS: %8.2f\n', Bias_OLS(d_idx)*100, Bias_FGLS(d_idx)*100);
        fprintf_both(fid, '       * MSE(beta)   | OLS: %8.2f | FGLS: %8.2f\n\n', MSE_OLS(d_idx)*100, MSE_FGLS(d_idx)*100);
    end
    
    res.d_grid = d_grid;
    res.Bias_OLS = Bias_OLS;
    res.MSE_OLS = MSE_OLS;
    res.Bias_FGLS = Bias_FGLS;
    res.MSE_FGLS = MSE_FGLS;
    res.MSE_S_OLS = MSE_S_OLS;
    res.MSE_S_FGLS = MSE_S_FGLS;
end

%% ========================================================================
%  EXPORT MODULE (Paper Exact Formatting)
% =========================================================================

function export_paper_exhibits_sec5_to_txt(All_Results, n_sims, folder_name, run_target, ts_file, ts_disp)
    filename = fullfile(folder_name, sprintf('Paper_Exhibits_SEC5_%dsims_%s.txt', n_sims, ts_file));
    fid = fopen(filename, 'w');
    if fid == -1, return; end
    
    fprintf(fid, '============================================================================================================\n');
    fprintf(fid, '               REPLICATION OF SECTION 5 PAPER EXHIBITS\n');
    fprintf(fid, '============================================================================================================\n');
    fprintf(fid, 'MC Iterations      : %d\n', n_sims);
    fprintf(fid, 'Date of Execution  : %s\n', ts_disp);
    fprintf(fid, 'Target Executed    : %s\n', upper(run_target));
    fprintf(fid, '============================================================================================================\n\n');
    
    run_tab67 = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab67');
    run_tab8  = strcmpi(run_target, 'all') || strcmpi(run_target, 'tab8');
    
    if run_tab67
        % --- EXPORT TABLE 6 ---
        fprintf(fid, '============================================================================================================\n');
        fprintf(fid, 'Table 6: MSE of \\hat{S}_{T,OLS} and \\hat{S}_{T,FGLS}\n');
        fprintf(fid, '============================================================================================================\n');
        write_tab67_block(fid, All_Results, 0.25, 0.6, true);
        write_tab67_block(fid, All_Results, 0.75, 0.6, true);
        write_tab67_block(fid, All_Results, 0.25, 0.8, true);
        write_tab67_block(fid, All_Results, 0.75, 0.8, true);
        fprintf(fid, '------------------------------------------------------------------------------------------------------------\n');
        fprintf(fid, 'Note: To facilitate readability each value is multiplied by 10^2.\n\n\n');
        
        % --- EXPORT TABLE 7 ---
        fprintf(fid, '============================================================================================================\n');
        fprintf(fid, 'Table 7: Bias and MSE of \\hat{\\beta}(\\hat{S}_{T,OLS}) and \\hat{\\beta}(\\hat{S}_{T,FGLS})\n');
        fprintf(fid, '============================================================================================================\n');
        write_tab67_block(fid, All_Results, 0.25, 0.6, false);
        write_tab67_block(fid, All_Results, 0.75, 0.6, false);
        write_tab67_block(fid, All_Results, 0.25, 0.8, false);
        write_tab67_block(fid, All_Results, 0.75, 0.8, false);
        fprintf(fid, '------------------------------------------------------------------------------------------------------------\n');
        fprintf(fid, 'Note: To facilitate readability each value is multiplied by 10^2.\n\n\n');
    end
    
    if run_tab8
        fprintf(fid, '============================================================================================================\n');
        fprintf(fid, 'Table 8: Bias and MSE of \\hat{\\beta}(\\hat{S}_{T,OLS}) and \\hat{\\beta}(\\hat{S}_{T,FGLS})\n');
        fprintf(fid, '============================================================================================================\n');
        write_tab8_block(fid, All_Results, 4);
        write_tab8_block(fid, All_Results, 8);
        fprintf(fid, '------------------------------------------------------------------------------------------------------------\n');
        fprintf(fid, 'Note: To facilitate readability each value is multiplied by 10^2.\n\n\n');
    end
    
    fclose(fid);
end

function write_tab67_block(fid, All_Results, rho, pi0, is_tab6)
    label_eu0 = sprintf('Tab67_rho%02d_pi%02d_eu00', round(rho*100), round(pi0*100));
    label_eu50 = sprintf('Tab67_rho%02d_pi%02d_eu50', round(rho*100), round(pi0*100));
    if ~isfield(All_Results, label_eu0) || ~isfield(All_Results, label_eu50), return; end
    res_0 = All_Results.(label_eu0);
    res_50 = All_Results.(label_eu50);
    
    fprintf(fid, '                                             rho_e = rho_u = 0                   rho_e = rho_u = 0.50\n');
    fprintf(fid, '                                     ---------------------------------   ---------------------------------\n');
    fprintf(fid, 'rho = %.2f, pi0 = %.1f                   d = 4       d = 16      d = 32       d = 4       d = 16      d = 32\n', rho, pi0);
    fprintf(fid, '------------------------------------------------------------------------------------------------------------\n');
    
    if is_tab6
        fprintf(fid, '%-36s', 'MSE(\hat{S}_{T,OLS})');
        fprintf(fid, '%12.2f%12.2f%12.2f', res_0.MSE_S_OLS * 100);
        fprintf(fid, '%12.2f%12.2f%12.2f\n', res_50.MSE_S_OLS * 100);
        
        fprintf(fid, '%-36s', 'MSE(\hat{S}_{T,FGLS})');
        fprintf(fid, '%12.2f%12.2f%12.2f', res_0.MSE_S_FGLS * 100);
        fprintf(fid, '%12.2f%12.2f%12.2f\n\n', res_50.MSE_S_FGLS * 100);
    else
        fprintf(fid, '%-36s', 'Bias(\hat{\beta}(\hat{S}_{T,OLS}))');
        fprintf(fid, '%12.2f%12.2f%12.2f', res_0.Bias_OLS * 100);
        fprintf(fid, '%12.2f%12.2f%12.2f\n', res_50.Bias_OLS * 100);
        
        fprintf(fid, '%-36s', 'Bias(\hat{\beta}(\hat{S}_{T,FGLS}))');
        fprintf(fid, '%12.2f%12.2f%12.2f', res_0.Bias_FGLS * 100);
        fprintf(fid, '%12.2f%12.2f%12.2f\n', res_50.Bias_FGLS * 100);
        
        fprintf(fid, '%-36s', 'MSE(\hat{\beta}(\hat{S}_{T,OLS}))');
        fprintf(fid, '%12.2f%12.2f%12.2f', res_0.MSE_OLS * 100);
        fprintf(fid, '%12.2f%12.2f%12.2f\n', res_50.MSE_OLS * 100);
        
        fprintf(fid, '%-36s', 'MSE(\hat{\beta}(\hat{S}_{T,FGLS}))');
        fprintf(fid, '%12.2f%12.2f%12.2f', res_0.MSE_FGLS * 100);
        fprintf(fid, '%12.2f%12.2f%12.2f\n\n', res_50.MSE_FGLS * 100);
    end
end

function write_tab8_block(fid, All_Results, d2_val)
    if d2_val == 4
        target = 'upper';
    else
        target = 'lower';
    end
    
    label_rho25 = sprintf('Tab8_%s_rho25', target);
    label_rho75 = sprintf('Tab8_%s_rho75', target);
    if ~isfield(All_Results, label_rho25) || ~isfield(All_Results, label_rho75), return; end
    res_25 = All_Results.(label_rho25);
    res_75 = All_Results.(label_rho75);
    
    fprintf(fid, '                                                rho = 0.25                          rho = 0.75\n');
    fprintf(fid, '                                     ---------------------------------   ---------------------------------\n');
    fprintf(fid, 'd_2 = %d, pi0 = 0.6                     d = 16       d = 24      d = 32      d = 16       d = 24      d = 32\n', d2_val);
    fprintf(fid, '------------------------------------------------------------------------------------------------------------\n');
    
    fprintf(fid, '%-36s', 'Bias(\hat{\beta}(\hat{S}_{T,OLS}))');
    fprintf(fid, '%12.2f%12.2f%12.2f', res_25.Bias_OLS * 100);
    fprintf(fid, '%12.2f%12.2f%12.2f\n', res_75.Bias_OLS * 100);
    
    fprintf(fid, '%-36s', 'Bias(\hat{\beta}(\hat{S}_{T,FGLS}))');
    fprintf(fid, '%12.2f%12.2f%12.2f', res_25.Bias_FGLS * 100);
    fprintf(fid, '%12.2f%12.2f%12.2f\n', res_75.Bias_FGLS * 100);
    
    fprintf(fid, '%-36s', 'MSE(\hat{\beta}(\hat{S}_{T,OLS}))');
    fprintf(fid, '%12.2f%12.2f%12.2f', res_25.MSE_OLS * 100);
    fprintf(fid, '%12.2f%12.2f%12.2f\n', res_75.MSE_OLS * 100);
    
    fprintf(fid, '%-36s', 'MSE(\hat{\beta}(\hat{S}_{T,FGLS}))');
    fprintf(fid, '%12.2f%12.2f%12.2f', res_25.MSE_FGLS * 100);
    fprintf(fid, '%12.2f%12.2f%12.2f\n\n', res_75.MSE_FGLS * 100);
end

%% ========================================================================
%  DGP FUNCTIONS (Verbatim Legacy Copied)
% =========================================================================

function [Y, D, Z, X, Index_S0] = generate_data_tab67(T, pi0, beta_dgp, d1, rho_v, rho_u, rho_e, burn_in, X_active, gamma1, gamma2, d_idx, beta_idx, i_sim)
    rng_seed = d_idx*1e6 + beta_idx*1e5 + i_sim;
    rng(rng_seed, "twister");
    
    theta_13 = d1 / sqrt(T);
    theta_2  = 0; 
    
    mu_vv = [0 0];
    Sigma_vv = [1, rho_v; rho_v, 1];
    
    Z = 1 + randn(T, 1);
    X = ones(T, 1); 
    D = zeros(T, 1); Y = zeros(T, 1);
    
    e_lag = 0; u_lag = 0;
    
    if rho_e ~= 0 || rho_u ~= 0
        for t_b = 1:burn_in
            innov = mvnrnd(mu_vv, Sigma_vv, 1);
            e_lag = rho_e * e_lag + innov(1);
            u_lag = rho_u * u_lag + innov(2);
        end
    end
    
    idx_r1_end = floor(T / 4);
    idx_r2_end = idx_r1_end + floor((1 - pi0) * T);
    
    for t = 1:T
        innov = mvnrnd(mu_vv, Sigma_vv, 1);
        e_t = rho_e * e_lag + innov(1);
        u_t = rho_u * u_lag + innov(2);
        
        if t <= idx_r1_end
            th = theta_13;
        elseif t <= idx_r2_end
            th = theta_2;
        else
            th = theta_13;
        end
        
        D(t) = th * Z(t) + gamma2 * X(t) + e_t;
        Y(t) = beta_dgp * D(t) + gamma1 * X(t) + u_t;
        
        e_lag = e_t; u_lag = u_t;
    end
    
    Index_S0 = ones(T, 1);
    Index_S0(idx_r1_end + 1 : idx_r2_end) = 0;
end

function [Y, D, Z, X, Index_S0] = generate_data_tab8(T, pi0, beta_dgp, d1, d2, rho_v, X_active, gamma1, gamma2, d_idx, beta_idx, i_sim)
    rng_seed = d_idx*1e6 + beta_idx*1e5 + i_sim;
    rng(rng_seed, "twister");
    
    theta_13 = d1 / sqrt(T);
    theta_2  = d2 / sqrt(T); 
    
    mu_eu = [0 0];
    Sigma_eu = [1, rho_v; rho_v, 1];
    
    Z = 1 + randn(T, 1);
    
    X = [];
    if X_active
        X = ones(T, 1); 
    end
    
    D = zeros(T, 1); Y = zeros(T, 1);
    
    idx_r1_end = floor(T / 4);
    idx_r2_end = idx_r1_end + floor((1 - pi0) * T);
    
    for t = 1:T
        R_err = mvnrnd(mu_eu, Sigma_eu, 1);
        e_D_t = R_err(1);
        u_Y_t = R_err(2);
        
        x_term_d = 0; x_term_y = 0;
        if X_active
            x_term_d = gamma2 * X(t,1);
            x_term_y = gamma1 * X(t,1);
        end
        
        if t <= idx_r1_end
            th = theta_13;
        elseif t <= idx_r2_end
            th = theta_2;
        else
            th = theta_13;
        end
        
        D(t) = th * Z(t) + x_term_d + e_D_t;
        Y(t) = beta_dgp * D(t) + x_term_y + u_Y_t;
    end
    
    Index_S0 = ones(T, 1);
    Index_S0(idx_r1_end + 1 : idx_r2_end) = 0;
end