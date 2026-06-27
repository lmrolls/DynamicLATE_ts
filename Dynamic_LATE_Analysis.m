function Results_Table = Dynamic_LATE_Analysis(Y, D, Z, X, options)
% DYNAMIC_LATE_ANALYSIS Performs structural extraction over a pi_L grid.
% Returns an iterative text output for the 3 phases, followed 
% by an asymptotic summary matrix conforming to econometric standards.
%
% INPUTS:
%   Y, D, Z, X : Data vectors and matrices (X can be [])
%   options    : Struct with piLow, m_breaks, alpha_level, beta_H0, etc.
%
% OUTPUT:
%   Results_Table : Struct containing the extracted fields for each pi_L.

    clc;
    warning('off','all');

    % Save original path for safe restoration
    original_path = path;
    cleanupObj = onCleanup(@() path(original_path));

    % --- 0. Setup and Parsing ---
    if nargin < 4; X = []; end
    if nargin < 5 || isempty(options); options = struct(); end
    
    if ~isfield(options, 'piLow'),       options.piLow = 0.5; end
    if ~isfield(options, 'm_breaks'),    options.m_breaks = 5; end
    if ~isfield(options, 'alpha_level'), options.alpha_level = 0.05; end
    if ~isfield(options, 'beta_H0'),     options.beta_H0 = 0.0; end
    if ~isfield(options, 'Hac_option'),  options.Hac_option = false; end
    
    T = size(Y, 1);
    q = size(Z, 2);
    
    % Generation of the sample fraction grid
    grid_piLow = round((options.piLow : 0.1 : 1.0) * 10) / 10;
    num_grid = length(grid_piLow);

    % Pre-allocation of vectors for the final matrix
    row_pi         = cell(1, num_grid);
    row_F_stat     = cell(1, num_grid);
    row_F_cv       = cell(1, num_grid);
    row_2sls_beta  = cell(1, num_grid);
    row_2sls_se    = cell(1, num_grid);
    row_2sls_t     = cell(1, num_grid);
    row_m_hat      = cell(1, num_grid);
    row_pi_hat     = cell(1, num_grid);
    row_dates      = cell(1, num_grid);
    
    row_AR_stat    = cell(1, num_grid);
    row_AR_pval    = cell(1, num_grid);
    row_LM_stat    = cell(1, num_grid);
    row_LM_pval    = cell(1, num_grid);
    row_LR_stat    = cell(1, num_grid);
    row_LR_pval    = cell(1, num_grid);

    % Add dependency paths
    addpath(fullfile(pwd, 'm_break'));
    addpath(fullfile(pwd, 'Section4', 'ftest'));
    addpath(fullfile(pwd, 'Section5'));
    addpath(fullfile(pwd, 'Section6'));

    fprintf('\n[+] Starting asymptotic extraction on pi_L grid (H0: beta = %.2f)...\n', options.beta_H0);

    %% CORE LOOP: SEQUENTIAL EVALUATION
    for i = 1:num_grid
        pi_val = grid_piLow(i);
        row_pi{i} = sprintf('pi_L = %.2f', pi_val);
        
        fprintf('\n===================================================================\n');
        fprintf(' ANALYSIS FOR pi_L = %.2f\n', pi_val);
        fprintf('===================================================================\n');
        
        % =================================================================
        % PHASE 1: First Stage Identification (F-Stat Optimization)
        % =================================================================
        test_options = struct('piLow', pi_val, 'm_breaks', options.m_breaks, 'Hac_option', options.Hac_option);
        F_Res = piLATE_Ftest(D, Z, X, test_options);
        
        cv_F = get_F_cv_Dynamic(pi_val, q, options.alpha_level);
        
        if F_Res.F_star > cv_F
            rej_str = 'Rej';
        else
            rej_str = 'Acc';
        end
        
        row_F_stat{i} = sprintf('%.4f', F_Res.F_star);
        row_F_cv{i}   = sprintf('[CV: %.2f] (%s)', cv_F, rej_str);
        row_m_hat{i}  = sprintf('%d', F_Res.best_m);
        row_pi_hat{i} = sprintf('%.2f', F_Res.best_pi);

        fprintf('>>> PHASE 1: F-STAT OPTIMIZATION\n');
        fprintf('    sup F* Stat       : %.4f\n', F_Res.F_star);
        fprintf('    CV (a=%.2f)       : %.2f (%s)\n', options.alpha_level, cv_F, rej_str);

        % =================================================================
        % PHASE 2: 2SLS Structural Estimation (Conditional on Phase 1)
        % =================================================================
        fprintf('\n>>> PHASE 2: 2SLS ESTIMATION\n');
        if ~isnan(F_Res.best_pi)
            [b_2sls, ~, se_2sls, t_2sls, dates_phase2] = est_2sls_pilate(Y, D, Z, X, F_Res.best_pi, F_Res.best_m);
            row_2sls_beta{i} = sprintf('%.4f', b_2sls);
            row_2sls_se{i}   = sprintf('(%.4f)', se_2sls);
            row_2sls_t{i}    = sprintf('[%.2f]', t_2sls);
            
            str_dates = sprintf('[%d,%d],', dates_phase2'); 
            if ~isempty(str_dates)
                str_dates(end) = []; 
            else
                str_dates = '[]';
            end
            
            row_dates{i} = str_dates;

            fprintf('    beta_hat          : %.4f\n', b_2sls);
            fprintf('    Standard Error    : %.4f\n', se_2sls);
            fprintf('    t-statistic       : %.4f\n', t_2sls);
            fprintf('    # of breaks (m_hat) : %d\n', F_Res.best_m);
            fprintf('    Fraction (pi_hat) : %.2f\n', F_Res.best_pi);
            fprintf('    Strongly-identified sub-sample dates (S_hat): %s\n', str_dates);
        else
            row_2sls_beta{i} = 'NaN'; 
            row_2sls_se{i}   = '(NaN)';
            row_2sls_t{i}    = '[NaN]';
            row_dates{i}     = '[NaN]';
            fprintf('    [!] Identification condition violated. Invalid sub-sample.\n');
        end

        % =================================================================
        % PHASE 3: Robust Inference (M2 Optimization)
        % =================================================================
        fprintf('\n>>> PHASE 3: IDENTIFICATION-ROBUST INFERENCE (H0: beta = %.2f)\n', options.beta_H0);
        params_core = struct('q', q, 'max_breaks', options.m_breaks, 'eps1', 0.05, ...
                             'vector_piL', (pi_val:0.1:1), 'step_pi', 0.05, ...
                             'Hac_option', options.Hac_option);
                             
        [~, stats_ST, cv_LR_struct] = robust_LATE_inference(Y, D, Z, X, options.beta_H0, params_core, []);
        
        cv_q = chi2inv(1 - options.alpha_level, q);
        cv_1 = chi2inv(1 - options.alpha_level, 1);
        
        if isfield(stats_ST, 'LM') && ~isnan(stats_ST.LM)
            AR_v = max(0, real(stats_ST.AR)); 
            p_AR = 1 - chi2cdf(AR_v * q, q);
            cv_AR_val = cv_q / q;
            
            LM_v = max(0, real(stats_ST.LM)); 
            p_LM = 1 - chi2cdf(LM_v, 1);
            cv_LM_val = cv_1;
            
            LR_v = max(0, real(stats_ST.LR)); 
            
            if isfield(cv_LR_struct, 'cv_ST') && ~isnan(cv_LR_struct.cv_ST)
                cv_LR_val = cv_LR_struct.cv_ST;
                p_LR = cv_LR_struct.pval_ST; 
            else
                cv_LR_val = cv_1;
                p_LR = NaN;
            end
            
            row_AR_stat{i} = sprintf('%.4f', AR_v);
            row_AR_pval{i} = sprintf('[p: %.3f] {CV: %.2f}', p_AR, cv_AR_val);
            row_LM_stat{i} = sprintf('%.4f', LM_v);
            row_LM_pval{i} = sprintf('[p: %.3f] {CV: %.2f}', p_LM, cv_LM_val);
            row_LR_stat{i} = sprintf('%.4f', LR_v);
            row_LR_pval{i} = sprintf('[p: %.3f] {CV: %.2f}', p_LR, cv_LR_val);
            
            fprintf('    AR Stat           : %.4f\n', AR_v);
            fprintf('    AR [pval] {CV}    : [p: %.3f] {CV: %.2f}\n', p_AR, cv_AR_val);
            fprintf('    LM Stat           : %.4f\n', LM_v);
            fprintf('    LM [pval] {CV}    : [p: %.3f] {CV: %.2f}\n', p_LM, cv_LM_val);
            fprintf('    LR Stat           : %.4f\n', LR_v);
            fprintf('    LR [pval] {CV}    : [p: %.3f] {CV: %.2f}\n', p_LR, cv_LR_val);
        else
            row_AR_stat{i}='NaN'; row_AR_pval{i}='[NaN] {NaN}';
            row_LM_stat{i}='NaN'; row_LM_pval{i}='[NaN] {NaN}';
            row_LR_stat{i}='NaN'; row_LR_pval{i}='[NaN] {NaN}';
            fprintf('    [!] M2 optimization failed. Insufficient rank or singular covariance.\n');
        end
    end

    %% CONSTRUCTION AND PRINT OF THE FINAL MATRIX
    cw = 26; 
    fmt = ['%-24s', repmat([' | %-', num2str(cw), 's'], 1, num_grid), '\n'];
    sep = repmat('-', 1, 24 + num_grid*(cw + 3));

    fprintf('\n\n===================================================================================================================\n');
    fprintf(' SUMMARY MATRIX\n');
    fprintf('===================================================================================================================\n');
    fprintf(fmt, 'METRIC \\ Grid', row_pi{:});
    fprintf('%s\n', sep);
    fprintf(fmt, 'sup F* Stat', row_F_stat{:});
    fprintf(fmt, '[CV] (Decision)', row_F_cv{:});
    fprintf(fmt, '2SLS \\beta', row_2sls_beta{:});
    fprintf(fmt, '(Std. Error)', row_2sls_se{:});
    fprintf(fmt, '[t-stat]', row_2sls_t{:});
    fprintf(fmt, '# of breaks (m_hat)', row_m_hat{:});
    fprintf(fmt, 'Fraction (pi_hat)', row_pi_hat{:});
    fprintf(fmt, 'Dates (S_hat)', row_dates{:});
    fprintf('%s\n', sep);
    fprintf(fmt, 'AR Stat', row_AR_stat{:});
    fprintf(fmt, '[p-value] {CV}', row_AR_pval{:});
    fprintf(fmt, 'LM Stat', row_LM_stat{:});
    fprintf(fmt, '[p-value] {CV}', row_LM_pval{:});
    fprintf(fmt, 'LR Stat', row_LR_stat{:});
    fprintf(fmt, '[p-value] {CV}', row_LR_pval{:});
    fprintf('%s\n\n', sep);

    Results_Table = struct();
    Results_Table.pi_Grid       = grid_piLow;
    Results_Table.F_Stat        = row_F_stat;
    Results_Table.F_CV          = row_F_cv;
    Results_Table.Beta_2SLS     = row_2sls_beta;
    Results_Table.SE_2SLS       = row_2sls_se;
    Results_Table.T_2SLS        = row_2sls_t;
    Results_Table.m_hat         = row_m_hat;
    Results_Table.pi_hat        = row_pi_hat;
    Results_Table.Dates         = row_dates;
    Results_Table.AR_Stat       = row_AR_stat;
    Results_Table.AR_pval       = row_AR_pval;
    Results_Table.LM_Stat       = row_LM_stat;
    Results_Table.LM_pval       = row_LM_pval;
    Results_Table.LR_Stat       = row_LR_stat;
    Results_Table.LR_pval       = row_LR_pval;

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function cv = get_F_cv_Dynamic(pi_L, q, alpha)
    % Extracts exact critical values (Table 4) for identification testing
    pi_grid = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    q_grid = [1, 2, 3, 4, 5, 10];
    
    pi_L = round(pi_L * 10) / 10;
    
    if alpha == 0.10
        table = [7.44, 5.18, 4.20, 3.70, 3.35, 2.53;
                 6.92, 4.76, 3.94, 3.46, 3.14, 2.40;
                 6.19, 4.44, 3.71, 3.28, 2.96, 2.31;
                 5.51, 4.02, 3.37, 2.79, 2.77, 2.18;
                 4.81, 3.59, 3.08, 2.78, 2.52, 2.04;
                 2.70, 2.32, 2.09, 1.94, 1.80, 1.59];
    elseif alpha == 0.01
        table = [12.27, 7.91, 6.08, 5.12, 4.56, 3.19;
                 11.63, 7.28, 5.73, 4.81, 4.37, 3.08;
                 10.94, 6.97, 5.56, 4.67, 4.19, 3.04;
                  9.73, 6.41, 5.06, 4.34, 3.89, 2.84;
                  8.68, 5.94, 4.62, 4.15, 3.65, 2.69;
                  6.68, 4.60, 3.70, 3.31, 2.99, 2.31];
    else % Default alpha = 0.05
        table = [8.90, 6.03, 4.75, 4.14, 3.74, 2.74;
                 8.28, 5.60, 4.49, 3.91, 3.51, 2.62;
                 7.55, 5.21, 4.26, 3.70, 3.31, 2.53;
                 6.84, 4.71, 3.83, 3.43, 3.13, 2.39;
                 6.04, 4.31, 3.58, 3.19, 2.89, 2.25;
                 3.85, 3.00, 2.57, 2.37, 2.16, 1.82];
    end
    
    pi_idx = find(abs(pi_grid - pi_L) < 1e-4, 1);
    q_idx = find(q_grid == q, 1);
    
    if isempty(pi_idx), pi_idx = 1; end % Fallback to 0.5
    if isempty(q_idx), q_idx = 1; end   % Fallback to q=1
    
    cv = table(pi_idx, q_idx);
end