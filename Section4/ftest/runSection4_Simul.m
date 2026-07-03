function Paper_Results = runSection4_Simul(n_sims_input, cv_reps_input, cv_option, run_target)
% RUNSECTION4_SIMUL Executes the simulation suite for exhibit replication (Section 4).
%
% SYNTAX:
%   res = runSection4_Simul()                         % Default: All (5000 sims, Fixed CV)
%   res = runSection4_Simul(5000, 10000, 2, 'tab5')   % Runs EXCLUSIVELY Table 5
%
% INPUTS:
%   n_sims_input  : Number of Monte Carlo simulations (Default: 5000)
%   cv_reps_input : Number of replications for the asymptotic CV (Default: 10000)
%   cv_option     : 1 = Asymptotically Simulated (via Compute_CVs)
%                   2 = Fixed Legacy (12.742 for F*, 3.84 for Full-Sample) - EXACT REPRODUCIBILITY
%   run_target    : Target exhibit to execute -> 'all' (default), 'fig5', 'tab5', 'fig6'
%
% =========================================================================

    clc; % DO NOT use 'clear', it destroys the function inputs
    restoredefaultpath;
    addpath('..\..\m_break'); 
    addpath('..\cvals'); 
    warning('off','all');

    % --- OPTIONAL INPUTS AND DEFAULT HANDLING ---
    if nargin < 1 || isempty(n_sims_input), n_sims = 5000; else, n_sims = n_sims_input; end
    if nargin < 2 || isempty(cv_reps_input), cv_reps = 10000; else, cv_reps = cv_reps_input; end
    if nargin < 3 || isempty(cv_option), cv_option = 2; end
    if nargin < 4 || isempty(run_target), run_target = 'all'; end
    
    valid_targets = {'all', 'fig5', 'tab5', 'fig6'};
    if ~ismember(lower(run_target), valid_targets)
        error('Error: run_target must be ''all'', ''fig5'', ''tab5'' or ''fig6''.');
    end

    piLow_param = 0.6; 

    % --- SETUP DESCRIPTIVE STRINGS FOR LOGS AND DIRECTORIES ---
    if cv_option == 1
        cv_str = 'SimCV';
        cv_desc = 'Simulated Asymptotic CV (via Compute_CVs)';
    else
        cv_str = 'FixedCV';
        cv_desc = 'Fixed Legacy CV (12.742 per pi-LATE F*, 3.84 per Full-Sample F)';
    end

    folder_name = sprintf('Paper_Exhibits_Results_%dsims_%s_%s', n_sims, cv_str, upper(run_target));
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    
    % --- TIMESTAMPS FOR SYNCHRONIZED LOGGING ---
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    timestamp_disp = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    log_filename = fullfile(folder_name, sprintf('Detailed_Logs_%dsims_%s_%s.txt', n_sims, cv_str, timestamp_file));
    fid_log = fopen(log_filename, 'w');
    if fid_log == -1, error('Unable to create the log file.'); end

    fprintf_both(fid_log, '===================================================================\n');
    fprintf_both(fid_log, ' STARTING SIMULATION SUITE: DYNAMIC LATE (SECTION 4)\n');
    fprintf_both(fid_log, ' Execution Target   : %s\n', upper(run_target));
    fprintf_both(fid_log, ' Execution Date/Time: %s\n', timestamp_disp);
    fprintf_both(fid_log, ' Critical Value Setup: %s\n', cv_desc);
    fprintf_both(fid_log, '===================================================================\n\n');

    fprintf_both(fid_log, 'Loading Empirical Data (Calibration)...\n');
    dta = readtable('dataForMatlabAltEstTick.csv');
    emp_data.FOMC = dta.FOMCused(401:800);
    emp_data.T_C = sum(emp_data.FOMC == 0);
    emp_data.T_P = sum(emp_data.FOMC == 1);
    emp_data.var_DNY2 = var(dta.DNY2(401:800));

    Paper_Results = struct();

    % --- STAGE 0: THEORETICAL CRITICAL VALUE INITIALIZATION ---
    fprintf_both(fid_log, '\n===================================================================\n');
    fprintf_both(fid_log, ' PHASE 0: Theoretical Critical Value Initialization\n');
    fprintf_both(fid_log, '===================================================================\n');
    if cv_option == 1
        fprintf_both(fid_log, '-> Computing Dynamic Asymptotic CV (piLow=%.2f, Reps=%d)...\n', piLow_param, cv_reps);
        cv_table = Compute_CVs(1, piLow_param, 1000, cv_reps, [], true);
        cv_teorico = cv_table.CV_95;
    else
        fprintf_both(fid_log, '-> Using Fixed Legacy CV...\n');
        cv_teorico = 12.742;
    end
    fprintf_both(fid_log, '   >> Assigned Theoretical CV F*: %.4f\n\n', cv_teorico);

    % =========================================================================
    % EXHIBIT 1: REPLICATION OF FIGURE 5 (Raw Power)
    % =========================================================================
    if strcmpi(run_target, 'all') || strcmpi(run_target, 'fig5')
        reset_parallel_pool(fid_log);
        
        fprintf_both(fid_log, '\n===================================================================\n');
        fprintf_both(fid_log, ' STARTING REPLICATION: FIGURE 5 (Raw Power, T=400)\n');
        fprintf_both(fid_log, '===================================================================\n');
        ar_grid_fig5 = [0, 0.25, 0.50, 0.75];
        for arCoef = ar_grid_fig5
            run_label = sprintf('Fig5_AR_%02d', arCoef*100);
            fprintf_both(fid_log, '\n>> Starting Iteration: %s\n', run_label);
            res = Run_MC_Scenario('Sec4', 0.6, arCoef, 0, n_sims, piLow_param, emp_data, cv_teorico);
            Paper_Results.(run_label) = res;
            print_scenario_results(fid_log, res, run_label);
        end
    end

    % =========================================================================
    % EXHIBIT 2: REPLICATION OF TABLE 5 (Size-Adjusted Power, Fixed Anomaly)
    % =========================================================================
    if strcmpi(run_target, 'all') || strcmpi(run_target, 'tab5')
        reset_parallel_pool(fid_log);
        
        fprintf_both(fid_log, '\n===================================================================\n');
        fprintf_both(fid_log, ' STARTING REPLICATION: TABLE 5 (Size-Adjusted Power, theta2 = -0.5)\n');
        fprintf_both(fid_log, '===================================================================\n');
        ar_grid_tab5 = [0.25, 0.75];
        for pi0 = [0.6, 0.8]
            for arCoef = ar_grid_tab5
                run_label = sprintf('Tab5_pi0_%02d_AR_%02d', pi0*10, arCoef*100);
                fprintf_both(fid_log, '\n>> Starting Iteration: %s\n', run_label);
                res = Run_MC_Scenario('M2_fixed', pi0, arCoef, 1, n_sims, piLow_param, emp_data, cv_teorico);
                Paper_Results.(run_label) = res;
                print_scenario_results(fid_log, res, run_label);
            end
        end
    end

    % =========================================================================
    % EXHIBIT 3: REPLICATION OF FIGURE 6 (Size-Adjusted Power, Proportional Anomaly)
    % =========================================================================
    if strcmpi(run_target, 'all') || strcmpi(run_target, 'fig6')
        reset_parallel_pool(fid_log);
        
        fprintf_both(fid_log, '\n===================================================================\n');
        fprintf_both(fid_log, ' STARTING REPLICATION: FIGURE 6 (Size-Adjusted Power, Proportional Anomaly)\n');
        fprintf_both(fid_log, '===================================================================\n');
        ar_grid_tab5 = [0.25, 0.75]; 
        for pi0 = [0.6, 0.8]
            for arCoef = ar_grid_tab5
                run_label = sprintf('Fig6_pi0_%02d_AR_%02d', pi0*10, arCoef*100);
                fprintf_both(fid_log, '\n>> Starting Iteration: %s\n', run_label);
                res = Run_MC_Scenario('M1_reverse', pi0, arCoef, 1, n_sims, piLow_param, emp_data, cv_teorico);
                Paper_Results.(run_label) = res;
                print_scenario_results(fid_log, res, run_label);
            end
        end
    end

    % --- TERMINAL SHUTDOWN ---
    reset_parallel_pool(fid_log); 
    fclose(fid_log);
    
    % Saving dated workspace file (.MAT)
    out_filename_mat = fullfile(folder_name, sprintf('Paper_Exhibits_Suite_%dsims_%s_%s.mat', n_sims, upper(run_target), timestamp_file));
    save(out_filename_mat, 'Paper_Results');
    
    % Exporting formatted text report (.TXT)
    export_paper_exhibits_to_txt(Paper_Results, n_sims, folder_name, cv_desc, timestamp_file, timestamp_disp);
    
    fprintf('\nExecution Completed Successfully (Target: %s). Outputs saved in: %s\\\n', upper(run_target), folder_name);
end 

%% ========================================================================
%  LOCAL FUNCTIONS (ORCHESTRATOR WRAPPER)
% =========================================================================

function reset_parallel_pool(fid)
    fprintf_both(fid, '\n[SYSTEM] Running memory cleanup and Parallel Pool reset...\n');
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
        fprintf_both(fid, '[SYSTEM] Previous Parallel Pool successfully closed.\n');
    else
        fprintf_both(fid, '[SYSTEM] No active Parallel Pool detected.\n');
    end
end

function Results = Run_MC_Scenario(dgp_type, pi0, arCoef, size_adjusted_flag, n_sims, piLow, emp_data, cv_teorico_piLATE)
    
    cv_teorico_FullSample = 3.84; 

    % --- EXACT REPLICATION OF SEEDS (arK MAPPING) ---
    if strcmp(dgp_type, 'Sec4')
        master_seed = 20250715; 
        T = 400; 
        thetaVecPower = 0.0:0.005:0.05;
        arVec_ref = [0, 0.25, 0.50, 0.75];
        arK = find(arVec_ref == arCoef);
    else
        master_seed = 20250808; 
        T = 200; 
        thetaVecPower = 0:2:24; 
        arK = 1; 
    end
    num_theta = length(thetaVecPower);

    Hac_option = (arCoef ~= 0); 

    % Settings configuration for the Core Engine API
    core_options = struct();
    core_options.piLow = piLow;
    core_options.Hac_option = Hac_option;
    core_options.m_breaks = 5;

    % --- STAGE 1: CRITICAL VALUE EVALUATION UNDER H0 (d = 0) ---
    FstatVec_piLATE_size = zeros(n_sims, 1);
    FstatVec_FullSample_size = zeros(n_sims, 1);

    parfor i_sims = 1:n_sims
        rng(master_seed + arK*1e6 + i_sims, "twister");
        
        [D, Z, ~, ~, X] = generate_Universal_DGP(T, pi0, 0, arCoef, dgp_type, 1, emp_data);
        
        % Core Engine Execution
        api_results = piLATE_Ftest(D, Z, X, core_options);
        
        FstatVec_piLATE_size(i_sims) = api_results.F_star;
        FstatVec_FullSample_size(i_sims) = api_results.F_FullSample;
    end

    Results.cv_piLATE = prctile(FstatVec_piLATE_size, 95);
    Results.cv_FullSample = prctile(FstatVec_FullSample_size, 95);
    
    Results.cv_teorico_piLATE = cv_teorico_piLATE;
    Results.cv_teorico_FullSample = cv_teorico_FullSample;
    Results.raw_size_piLATE = mean(FstatVec_piLATE_size > cv_teorico_piLATE);
    Results.raw_size_FullSample = mean(FstatVec_FullSample_size > cv_teorico_FullSample);

    % --- STAGE 2: POWER CURVE COMPUTATION ---
    if size_adjusted_flag
        cv_final_piLATE = Results.cv_piLATE;
        cv_final_FullSample = Results.cv_FullSample;
    else
        cv_final_piLATE = cv_teorico_piLATE;
        cv_final_FullSample = cv_teorico_FullSample;
    end

    FstatVec_piLATE_power_Matrix = zeros(n_sims, num_theta);
    FstatVec_FullSample_power_Matrix = zeros(n_sims, num_theta);

    parfor i_sims = 1:n_sims
        temp_pL = zeros(1, num_theta);
        temp_FS = zeros(1, num_theta);
        
        for thetaK = 1:num_theta
            theta_val_power = thetaVecPower(thetaK);
            rng(master_seed + 1e8 + arK*1e6 + thetaK*1e4 + i_sims, "twister");
            
            [D, Z, ~, ~, X] = generate_Universal_DGP(T, pi0, theta_val_power, arCoef, dgp_type, 0, emp_data);
            
            % Core Engine Execution
            api_results = piLATE_Ftest(D, Z, X, core_options);
            
            temp_pL(thetaK) = api_results.F_star;
            temp_FS(thetaK) = api_results.F_FullSample;
        end
        FstatVec_piLATE_power_Matrix(i_sims, :) = temp_pL;
        FstatVec_FullSample_power_Matrix(i_sims, :) = temp_FS;
    end
    
    FstatTotal_Power_piLATE = zeros(1, num_theta);
    FstatTotal_Power_FullSample = zeros(1, num_theta);
    
    for thetaK = 1:num_theta
        FstatTotal_Power_piLATE(1, thetaK) = mean(FstatVec_piLATE_power_Matrix(:, thetaK) > cv_final_piLATE);
        FstatTotal_Power_FullSample(1, thetaK) = mean(FstatVec_FullSample_power_Matrix(:, thetaK) > cv_final_FullSample);
    end

    Results.power_piLATE = FstatTotal_Power_piLATE;
    Results.power_FullSample = FstatTotal_Power_FullSample;
    Results.thetaVecPower = thetaVecPower;
    Results.dgp_type = dgp_type;
    Results.arCoef = arCoef;
    Results.pi0 = pi0;
end

function y = AR1_Sigma(T, alpha, sigma)
    raw_e = randn(T*2, 1);
    raw_y0 = randn();
    e = sigma * raw_e;
    y = zeros(T*2, 1);
    y(1,1) = sigma * raw_y0; 
    for i = 2:T*2, y(i,1) = alpha * y(i-1,1) + e(i,1); end
    y = y(T+1:end, 1);
end

function [D, Z, q, p, X] = generate_Universal_DGP(T, pi0, d_val, arCoef, dgp_type, is_H0, emp_data)
    if strcmp(dgp_type, 'Sec4')
        theta_dgp = d_val; 
        Z_base = zeros(T, 1);
        Z_base(emp_data.FOMC == 1) = T / emp_data.T_P;
        Z_base(emp_data.FOMC == 0) = -T / emp_data.T_C;
        q = 1; p = 0; X = [];
        sigmaE = sqrt(emp_data.var_DNY2);
        
        if arCoef == 0, e_D = sigmaE * randn(T, 1); else, e_D = AR1_Sigma(T, arCoef, sigmaE); end
        
        D = zeros(T, 1);
        idx_r1 = floor(T/4);
        idx_r2 = idx_r1 + floor((1 - pi0) * T);
        
        for t = 1:T
            if t <= idx_r1
                D(t) = e_D(t) / (1 - theta_dgp * Z_base(t));
            elseif t <= idx_r2
                D(t) = e_D(t) / (1 - 0 * Z_base(t));
            else
                D(t) = e_D(t) / (1 - theta_dgp * Z_base(t));
            end
        end
        Z = Z_base .* D;
    else
        theta_dgp = d_val / sqrt(T);
        Z_instrument = 1 + sqrt(1)*randn(T, 1); 
        q = 1; X = ones(T, 1); p = 1; gamma2 = 1; sigmaE = 1;
        
        if arCoef == 0, e_D = sigmaE * randn(T, 1); else, e_D = AR1_Sigma(T, arCoef, sigmaE); end
        
        D = zeros(T, 1);
        idx_r1 = floor(T / 4);
        idx_r2 = idx_r1 + floor((1 - pi0) * T);
        
        for t = 1:T
            x_term = gamma2 * X(t);
            if t <= idx_r1
                D(t) = theta_dgp * Z_instrument(t) + x_term + e_D(t);
            elseif t <= idx_r2
                if is_H0
                    D(t) = theta_dgp * Z_instrument(t) + x_term + e_D(t);
                else
                    if strcmp(dgp_type, 'M1_reverse')
                        D(t) = -(theta_dgp / 2) * Z_instrument(t) + x_term + e_D(t);
                    elseif strcmp(dgp_type, 'M2_fixed')
                        D(t) = -0.5 * Z_instrument(t) + x_term + e_D(t);
                    end
                end
            else
                D(t) = theta_dgp * Z_instrument(t) + x_term + e_D(t);
            end
        end
        Z = Z_instrument;
    end
end

function fprintf_both(fid, fmt, varargin)
    fprintf(1, fmt, varargin{:});
    if fid ~= -1, fprintf(fid, fmt, varargin{:}); end
end

function print_scenario_results(fid, res, run_name)
    fprintf_both(fid, '--- Analytical Results for Scenario [%s | arCoef=%.2f | pi0=%.1f] ---\n', res.dgp_type, res.arCoef, res.pi0);
    fprintf_both(fid, '  > Nominal size evaluated against CV=%.4f (F*) and CV=%.4f (Full) \n', res.cv_teorico_piLATE, res.cv_teorico_FullSample);
    fprintf_both(fid, '  > Raw Empirical Size -> F* : %.4f | Full-Sample : %.4f\n', res.raw_size_piLATE, res.raw_size_FullSample);
    fprintf_both(fid, '  > New Empirical CVs (95%%)-> F* : %.4f | Full-Sample : %.4f\n', res.cv_piLATE, res.cv_FullSample);
    
    if strcmp(res.dgp_type, 'Sec4')
        fprintf_both(fid, '\n  [ RAW POWER TABLE ]\n');
        param_label = 'Theta';
    else
        fprintf_both(fid, '\n  [ SIZE-ADJUSTED POWER TABLE ]\n');
        param_label = 'd';
    end
    
    fprintf_both(fid, '  --------------------------------------------------------\n');
    fprintf_both(fid, '  %10s | %15s | %15s\n', param_label, 'pi-LATE F*', 'Full-Sample F');
    fprintf_both(fid, '  --------------------------------------------------------\n');
    for k = 1:length(res.thetaVecPower)
        fprintf_both(fid, '  %10.4f | %15.4f | %15.4f\n', res.thetaVecPower(k), res.power_piLATE(k), res.power_FullSample(k));
    end
    fprintf_both(fid, '  --------------------------------------------------------\n\n');
end

function export_paper_exhibits_to_txt(All_Results, n_sims, folder_name, cv_desc, ts_file, ts_disp)
    filename = fullfile(folder_name, sprintf('Paper_Exhibits_%dsims_%s.txt', n_sims, ts_file));
    fid = fopen(filename, 'w');
    if fid == -1, return; end
    
    fprintf(fid, '====================================================================\n');
    fprintf(fid, '               PAPER EXHIBIT REPLICATION                            \n');
    fprintf(fid, '====================================================================\n');
    fprintf(fid, 'MC Iterations      : %d\n', n_sims);
    fprintf(fid, 'Date of Execution  : %s\n', ts_disp);
    fprintf(fid, 'Critical Val Logic : %s\n', cv_desc);
    fprintf(fid, '====================================================================\n\n');
    fields = fieldnames(All_Results);
    
    % --- TABLE 5 ---
    tab5_fields = fields(startsWith(fields, 'Tab5_'));
    if ~isempty(tab5_fields)
        fprintf(fid, '====================================================================\n');
        fprintf(fid, ' TABLE 5: Size and size-adjusted power with fixed anomaly\n');
        fprintf(fid, '====================================================================\n');
        for i = 1:length(tab5_fields)
            res = All_Results.(tab5_fields{i});
            d_vec = res.thetaVecPower;
            if res.pi0 == 0.6, d_target = [0, 4, 8, 12, 16, 20, 24]; else, d_target = [0, 2, 6, 10, 14, 18, 22]; end
            [~, print_indices] = ismember(d_target, d_vec);
            print_indices = print_indices(print_indices > 0);
            
            fprintf(fid, '>>> [rho = %.2f, pi0 = %.1f]\n', res.arCoef, res.pi0);
            fprintf(fid, '%-10s | %-15s | %-15s\n', 'd', 'pi-LATE F*', 'Full-Sample');
            fprintf(fid, '----------------------------------------------\n');
            for k = print_indices
                pL = res.power_piLATE(k); if d_vec(k) == 0, pL = res.raw_size_piLATE; end
                pFS = res.power_FullSample(k); if d_vec(k) == 0, pFS = res.raw_size_FullSample; end
                fprintf(fid, '%-10.4f | %-15.4f | %-15.4f\n', d_vec(k), pL, pFS);
            end
            fprintf(fid, '\n');
        end
    end
    
    % --- FIGURE 6 ---
    fig6_fields = fields(startsWith(fields, 'Fig6_'));
    if ~isempty(fig6_fields)
        fprintf(fid, '\n====================================================================\n');
        fprintf(fid, ' FIGURE 6: Size-Adjusted Power for proportional anomaly\n');
        fprintf(fid, '====================================================================\n');
        for i = 1:length(fig6_fields)
            res = All_Results.(fig6_fields{i});
            fprintf(fid, '>>> [rho = %.2f, pi0 = %.1f]\n', res.arCoef, res.pi0);
            fprintf(fid, '%-10s | %-15s | %-15s\n', 'd', 'pi-LATE F*', 'Full-Sample');
            fprintf(fid, '----------------------------------------------\n');
            print_indices = 1:2:length(res.thetaVecPower);
            for k = print_indices
                fprintf(fid, '%-10.4f | %-15.4f | %-15.4f\n', res.thetaVecPower(k), res.power_piLATE(k), res.power_FullSample(k));
            end
            fprintf(fid, '\n');
        end
    end

    % --- FIGURE 5 ---
    fig5_fields = fields(startsWith(fields, 'Fig5_'));
    if ~isempty(fig5_fields)
        fprintf(fid, '====================================================================\n');
        fprintf(fid, ' FIGURE 5: Power curves for calibrated DGP (Raw Power)\n');
        fprintf(fid, '====================================================================\n');
        for i = 1:length(fig5_fields)
            res = All_Results.(fig5_fields{i});
            fprintf(fid, '>>> [rho = %.2f, Calibration T=400]\n', res.arCoef);
            fprintf(fid, '%-10s | %-15s | %-15s\n', 'Theta', 'pi-LATE F*', 'Full-Sample');
            fprintf(fid, '----------------------------------------------\n');
            print_indices = 1:2:length(res.thetaVecPower);
            for k = print_indices
                fprintf(fid, '%-10.4f | %-15.4f | %-15.4f\n', res.thetaVecPower(k), res.power_piLATE(k), res.power_FullSample(k));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
end