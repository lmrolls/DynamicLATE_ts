function test_Dynamic_LATE_all_DGPs()
% TEST_DYNAMIC_LATE_ALL_DGPS
% Executes a sequential test of the Dynamic_LATE_Analysis pipeline on
% the three main Data Generating Processes (DGPs) of the paper:
% - Scenario 1 (Section 5): Exogenous covariates, i.i.d. errors.
% - Scenario 2 (Section 6): Perverse identification (local weak), i.i.d. errors.

    clc;
    close all;

    fprintf('===================================================================\n');
    fprintf(' UNIFIED PIPELINE TEST: DYNAMIC LATE ANALYSIS (SECTIONS 4, 5, 6) \n');
    fprintf('===================================================================\n');

    % Shared baseline options
    options_base = struct();
    options_base.piLow = 0.6;
    options_base.m_breaks = 5;
    options_base.alpha_level = 0.05;


    %% =====================================================================
    % SCENARIO 1: SECTION 5 (Structural Estimation with Exogenous X)
    % =====================================================================
    fprintf('\n>>> SCENARIO 1: SECTION 5 (Exogenous Covariates X, HAC = OFF)\n');
    T5 = 200; pi0_5 = 0.6; rho_eu_5 = 0.25; d_val_5 = 16; beta_true_5 = 1.0;
    
    options_5 = options_base;
    options_5.Hac_option = false; % i.i.d. errors, White-robust (HC1) is sufficient
    options_5.beta_H0 = 1.0;      % Centering H0 on the true beta to verify size
    
    disp('Press any key to proceed to Scenario 1...');
    pause;

    rng_seed_5 = 202605;
    [Y5, D5, Z5, X5] = dgp_sec5(T5, pi0_5, d_val_5, rho_eu_5, beta_true_5, rng_seed_5);
    
    fprintf(' Section 5 Data Generation completed. Starting analysis...\n');
    Results_Sec5 = Dynamic_LATE_Analysis(Y5, D5, Z5, X5, options_5);


    %% =====================================================================
    % SCENARIO 2: SECTION 6 (Perverse Identification / Local Weak)
    % =====================================================================
    fprintf('\n>>> SCENARIO 2: SECTION 6 (Perverse Identification, HAC = OFF)\n');
    T6 = 200; pi0_6 = 0.6; rho_eu_6 = 0.25; d_val_6 = 24; beta_true_6 = 0.0;
    
    options_6 = options_base;
    options_6.Hac_option = false; % i.i.d. errors, White-robust (HC1) is sufficient
    options_6.beta_H0 = 0.0;
    
    disp('Press any key to proceed to Scenario 2...');
    pause;

    rng_seed_6 = 202606;
    [Y6, D6, Z6, X6] = dgp_sec6(T6, pi0_6, d_val_6, rho_eu_6, beta_true_6, rng_seed_6);
    
    fprintf(' Section 6 Data Generation completed. Starting analysis...\n');
    Results_Sec6 = Dynamic_LATE_Analysis(Y6, D6, Z6, X6, options_6);
    
    fprintf('\n===================================================================\n');
    fprintf(' SEQUENTIAL TEST COMPLETED. Structures saved in the workspace.\n');
    fprintf('===================================================================\n');
end

%% ========================================================================
%  LOCAL DGP HELPERS
% =========================================================================

function [Y, D, Z, X] = dgp_sec5(T, pi0, d_val, rho_eu, beta_true, rng_seed)
    rng(rng_seed, "twister");
    theta_dgp = d_val / sqrt(T);
    gamma1 = 1; gamma2 = 1;
    
    Z = 1 + randn(T, 1);
    X = ones(T, 1);
    D = zeros(T, 1); 
    Y = zeros(T, 1);
    
    Sigma_eu = [1 rho_eu; rho_eu 1];
    idx_r1 = floor(T/4);
    idx_r2 = idx_r1 + floor((1 - pi0) * T);
    
    for t = 1:T
        err = mvnrnd([0 0], Sigma_eu, 1);
        if t <= idx_r1 || t > idx_r2
            D(t) = theta_dgp * Z(t) + gamma2 * X(t) + err(1);
        else
            D(t) = 0 * Z(t) + gamma2 * X(t) + err(1); % Fixed anomaly
        end
        Y(t) = beta_true * D(t) + gamma1 * X(t) + err(2);
    end
end

function [Y, D, Z, X] = dgp_sec6(T, pi0, d_val, rho_eu, beta_true, rng_seed)
    rng(rng_seed, "twister");
    theta_strong = d_val / sqrt(T);
    theta_weak   = -0.5 * theta_strong; % Perverse identification
    
    Z = 1 + randn(T, 1);
    X = []; % M1
    D = zeros(T, 1); 
    Y = zeros(T, 1);
    
    Sigma_eu = [1 rho_eu; rho_eu 1];
    idx_r1 = floor(T/4);
    idx_r2 = idx_r1 + floor((1 - pi0) * T);
    
    for t = 1:T
        err = mvnrnd([0 0], Sigma_eu, 1);
        if t <= idx_r1 || t > idx_r2
            D(t) = theta_strong * Z(t) + err(1);
        else
            D(t) = theta_weak * Z(t) + err(1);
        end
        Y(t) = beta_true * D(t) + err(2);
    end
end