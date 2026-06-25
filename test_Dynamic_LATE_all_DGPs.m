function test_Dynamic_LATE_all_DGPs()
% TEST_DYNAMIC_LATE_ALL_DGPS
% Executes a sequential test of the Dynamic_LATE_Analysis pipeline on
% the three main Data Generating Processes (DGPs) of the paper:
% - Scenario 1 (Section 4): AR(1) errors, discrete instrument.
% - Scenario 2 (Section 5): Exogenous covariates, i.i.d. errors.
% - Scenario 3 (Section 6): Perverse identification (local weak), i.i.d. errors.

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
    % SCENARIO 1: SECTION 4 (Empirical Calibration, AR(1) Errors)
    % =====================================================================
    fprintf('\n>>> SCENARIO 1: SECTION 4 (AR(1) Errors, HAC = ON)\n');
    T4 = 400; pi0_4 = 0.6; arCoef_4 = 0.25; d_val_4 = 0.02; beta_true_4 = 0.0;
    
    options_4 = options_base;
    options_4.Hac_option = true; % Mandatory to correct for AR(1) persistence
    options_4.beta_H0 = 0.0;
    
    rng_seed_4 = 202604;
    [Y4, D4, Z4, X4] = dgp_sec4(T4, pi0_4, d_val_4, arCoef_4, beta_true_4, rng_seed_4);
    
    fprintf(' Section 4 Data Generation completed. Starting analysis...\n');
    Results_Sec4 = Dynamic_LATE_Analysis(Y4, D4, Z4, X4, options_4);
    
    disp('Press any key to proceed to Scenario 2...');
    pause;

    %% =====================================================================
    % SCENARIO 2: SECTION 5 (Structural Estimation with Exogenous X)
    % =====================================================================
    fprintf('\n>>> SCENARIO 2: SECTION 5 (Exogenous Covariates X, HAC = OFF)\n');
    T5 = 200; pi0_5 = 0.6; rho_eu_5 = 0.25; d_val_5 = 16; beta_true_5 = 1.0;
    
    options_5 = options_base;
    options_5.Hac_option = false; % i.i.d. errors, White-robust (HC1) is sufficient
    options_5.beta_H0 = 1.0;      % Centering H0 on the true beta to verify size
    
    rng_seed_5 = 202605;
    [Y5, D5, Z5, X5] = dgp_sec5(T5, pi0_5, d_val_5, rho_eu_5, beta_true_5, rng_seed_5);
    
    fprintf(' Section 5 Data Generation completed. Starting analysis...\n');
    Results_Sec5 = Dynamic_LATE_Analysis(Y5, D5, Z5, X5, options_5);

    disp('Press any key to proceed to Scenario 3...');
    pause;

    %% =====================================================================
    % SCENARIO 3: SECTION 6 (Perverse Identification / Local Weak)
    % =====================================================================
    fprintf('\n>>> SCENARIO 3: SECTION 6 (Perverse Identification, HAC = OFF)\n');
    T6 = 200; pi0_6 = 0.6; rho_eu_6 = 0.25; d_val_6 = 24; beta_true_6 = 0.0;
    
    options_6 = options_base;
    options_6.Hac_option = false; % i.i.d. errors, White-robust (HC1) is sufficient
    options_6.beta_H0 = 0.0;
    
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

function [Y, D, Z, X] = dgp_sec4(T, pi0, d_val, arCoef, beta_true, rng_seed)
    rng(rng_seed, "twister");
    theta_dgp = d_val; 
    
    % Mock empirical calibration (simulates FOMC)
    T_P = floor(T/2); T_C = T - T_P;
    FOMC = double(rand(T, 1) > 0.5);
    Z_base = zeros(T, 1);
    Z_base(FOMC == 1) = T / T_P;
    Z_base(FOMC == 0) = -T / T_C;
    
    X = []; % M1: No exogenous variables
    sigmaE = 1.0; % Standardized unit variance
    
    % AR(1) perturbations generation
    raw_e = randn(T*2, 1);
    e_D = zeros(T*2, 1);
    e_D(1) = sigmaE * randn(); 
    for t = 2:T*2, e_D(t) = arCoef * e_D(t-1) + sigmaE * raw_e(t); end
    e_D = e_D(T+1:end);
    
    D = zeros(T, 1);
    idx_r1 = floor(T/4);
    idx_r2 = idx_r1 + floor((1 - pi0) * T);
    
    for t = 1:T
        if t <= idx_r1 || t > idx_r2
            D(t) = e_D(t) / (1 - theta_dgp * Z_base(t));
        else
            D(t) = e_D(t); % Structurally weak instrument (th=0)
        end
    end
    Z = Z_base .* D;
    Y = beta_true * D + sigmaE * randn(T, 1);
end

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