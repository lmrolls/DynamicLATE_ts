function [beta_hat, S_vec, SE_robust, t_stat, dates] = est_2sls_pilate(Y, D, Z, X, pi0, m, Hac_option)
% EST_2SLS_PILATE Computes the 2SLS estimator on the optimal sub-sample identified
% via the Bai-Perron procedure (Phase 1).
%
% INPUTS:
%   Y, D, Z, X : Original data vectors and matrices
%   pi0        : Sample fraction
%   m          : Number of structural breaks (e.g., 5)
%   Hac_option : (Optional) Boolean. If true, applies Newey-West HAC.
%                If false or omitted, applies White HC1.

    % --- BACKWARD COMPATIBILITY BLOCK ---
    % If invoked by Section 5 simulations (6 arguments), force HC1
    if nargin < 7 || isempty(Hac_option)
        Hac_option = false; 
    end
    T = size(Y, 1);
    q = size(Z, 2);
    p = size(X, 2);
    h = round(0.1 * T); 
    
    [~, datevec, ~] = Procedure_piL_Full_Sample(T, D, Z, q, m, h, 0.1, 1, 0, 1, 1, 1, 0, X, p, 1e-4, 20, 0, 0, 1, 1);
    dates = S_Dates_m(D, Z, X, datevec, T, pi0, m);
    
    [beta_hat, S_vec, SE_robust, t_stat] = internal_ols_calc(Y, D, Z, X, dates, T, Hac_option);
end

function [beta_val, S_vec, SE_val, t_stat_val] = internal_ols_calc(Y, D, Z, X, Sup_S, T, Hac_option)
    S_vec = zeros(T, 1);
    for irow=1:size(Sup_S,1), S_vec(Sup_S(irow,1):Sup_S(irow,2))=1; end
    
    idx = (S_vec == 1);
    Y_c = Y(idx);
    D_c = D(idx);
    Z_c = Z(idx, :);
    X_c = []; 
    if size(X,2) > 0, X_c = X(idx, :); end
    
    q = size(Z,2); p = size(X,2);
    T_eff = sum(S_vec); 
    
    if T_eff <= (q+p)
        beta_val = NaN; SE_val = NaN; t_stat_val = NaN; return; 
    end
    
    Instruments_c = [Z_c, X_c];
    P_Instruments_c = Instruments_c * (Instruments_c' * Instruments_c)^-1 * Instruments_c';
    D_hat_c = P_Instruments_c * D_c;
    
    % Rigorous algebraic separation: Projected Regressors vs Observed Regressors
    Regressors_proj = [D_hat_c, X_c];
    Regressors_obs  = [D_c, X_c];
    
    % Passing the HAC parameter to the covariance calculation subroutine
    [beta_hat_struct, SE_robust_struct, ~, ~, t_stat_struct] = OLS_Robust_2SLS(Y_c, Regressors_proj, Regressors_obs, T_eff, Hac_option);
    
    beta_val = beta_hat_struct(1);
    SE_val = SE_robust_struct(1);
    t_stat_val = t_stat_struct(1);
end

function [beta_hat, SE_robust, CI, resid, t_stat] = OLS_Robust_2SLS(Y, X_proj, X_obs, T_eff, Hac_option)
% Computes the 2SLS estimator and asymptotic inference (HAC or White-robust).
% Employs X_proj for point estimation and X_obs for structural residuals.
    [~, nvar] = size(X_proj);
    if size(Y, 2) > 1, Y = Y'; end
    % 1. Point Estimation (HAC-invariant)
    beta_hat = (X_proj' * X_proj) \ (X_proj' * Y);
    % 2. Corrected Structural Residuals (Based on observed endogenous variables)
    resid = Y - X_obs * beta_hat;
    % 3. Asymptotic Variance
    XX_inv = (X_proj' * X_proj) \ eye(nvar); 
    
    if Hac_option
        Lag = floor(T_eff^(1/3));
        % NeweyWest computes the sum of cross-products and already divides by (T_eff - nvar).
        % Hence J_hac = Q_sum / (T_eff - nvar).
        J_hac = NeweyWest(resid, X_proj, Lag, 0);
        
        % The sandwich HC1 formula is: (T_eff / (T_eff - nvar)) * XX_inv * Q_sum * XX_inv
        % Substituting Q_sum = J_hac * (T_eff - nvar), the degrees of freedom simplify.
        % The resulting scaling factor multiplying the structure becomes exactly T_eff.
        VCE_robust = T_eff * (XX_inv * J_hac * XX_inv);
    else
        % --- LEGACY BLOCK: Bit-exact reproduction for Section 5 (HC1 White) ---
        X_res = X_proj .* resid; 
        meat = X_res' * X_res;
        scale = T_eff / (T_eff - nvar);
        VCE_robust = scale * (XX_inv * meat * XX_inv);
    end
    % 4. Inference
    SE_robust = sqrt(diag(VCE_robust));
    t_stat = beta_hat ./ SE_robust;
    CI = [beta_hat - 1.96 * SE_robust, beta_hat + 1.96 * SE_robust];
end