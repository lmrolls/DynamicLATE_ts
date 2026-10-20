function [beta_hat, SE_robust, t_stat] = est_iv_fullsample(Y, D, Z, X, Hac_option)
% EST_IV_FULLSAMPLE Full-sample IV (2SLS) estimator of beta, i.e. the rival
% benchmark \hat{beta}_{FS} of Section 5 (Tables 7 and 8). No sub-sample
% selection is performed: all T observations are used.
%
% SYNTAX:
%   [beta_hat, SE_robust, t_stat] = est_iv_fullsample(Y, D, Z, X)
%   [beta_hat, SE_robust, t_stat] = est_iv_fullsample(Y, D, Z, X, Hac_option)
%
% INPUTS:
%   Y, D       : Outcome and treatment (T x 1)
%   Z, X       : Instruments (T x q) and exogenous regressors (T x p, may be [])
%   Hac_option : (Optional) Boolean. true  -> Newey-West HAC standard errors
%                (Lag = floor(T^(1/3)); requires NeweyWest.m from m_break on the path);
%                false or omitted -> White HC1 standard errors.
%
% OUTPUTS:
%   beta_hat, SE_robust, t_stat : coefficient on D (first element)
%
% NOTES:
%   * Point estimate: literal port of estimate_full_sample / OLS_Robust_full in the
%     original Section 5 replication package (same projection expression and
%     backslash solve, hence bit-exact with the published \hat{beta}_{FS}).
%   * Standard errors follow est_2sls_pilate.m: structural residuals computed on
%     the observed regressors [D, X], not on the projected ones.
%   * With identical regressors in the two equations, full-sample SUR/FGLS reduces
%     to equation-by-equation OLS (Zellner), so beta_hat is also the full-sample
%     counterpart of the FGLS estimator: a single rival benchmark covers both.
%   * Used only by the Section 5 simulation script; not part of the standalone
%     pi-LATE functions or of Dynamic_LATE_Analysis.m.

    if nargin < 5 || isempty(Hac_option)
        Hac_option = false;
    end

    T = size(Y, 1);
    q = size(Z, 2);
    p = size(X, 2);

    if T <= (q + p)
        beta_hat = NaN; SE_robust = NaN; t_stat = NaN;
        return;
    end

    % First stage: projection of D on the full-sample instrument set [Z, X]
    Instruments = [Z, X];
    P_Instruments = Instruments * (Instruments' * Instruments)^-1 * Instruments';
    D_hat = P_Instruments * D;

    % Second stage: Y on the projected regressors; residuals on the observed ones
    Regressors_proj = [D_hat, X];
    Regressors_obs  = [D, X];

    [beta_hat_struct, SE_robust_struct, t_stat_struct] = ols_robust_2sls_local(Y, Regressors_proj, Regressors_obs, T, Hac_option);

    beta_hat  = beta_hat_struct(1);
    SE_robust = SE_robust_struct(1);
    t_stat    = t_stat_struct(1);
end

function [beta_hat, SE_robust, t_stat] = ols_robust_2sls_local(Y, X_proj, X_obs, T_eff, Hac_option)
% 2SLS point estimate and robust (HC1 or HAC) inference. Same algebra as
% OLS_Robust_2SLS in est_2sls_pilate.m, replicated here so that the helper is
% self-contained (NeweyWest.m is needed only when Hac_option is true).
    [~, nvar] = size(X_proj);
    if size(Y, 2) > 1, Y = Y'; end

    % 1. Point estimation (identical to OLS_Robust_full: (X'X)\(X'Y))
    beta_hat = (X_proj' * X_proj) \ (X_proj' * Y);

    % 2. Structural residuals (observed endogenous regressors)
    resid = Y - X_obs * beta_hat;

    % 3. Asymptotic variance
    XX_inv = (X_proj' * X_proj) \ eye(nvar);
    if Hac_option
        Lag = floor(T_eff^(1/3));
        J_hac = NeweyWest(resid, X_proj, Lag, 0);     % = Q_sum / (T_eff - nvar)
        VCE_robust = T_eff * (XX_inv * J_hac * XX_inv);
    else
        X_res = X_proj .* resid;
        meat  = X_res' * X_res;
        scale = T_eff / (T_eff - nvar);
        VCE_robust = scale * (XX_inv * meat * XX_inv);
    end

    % 4. Inference
    SE_robust = sqrt(diag(VCE_robust));
    t_stat = beta_hat ./ SE_robust;
end
