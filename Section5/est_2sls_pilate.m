function [beta_hat, S_vec, SE_robust, t_stat, dates] = est_2sls_pilate(Y, D, Z, X, pi0, m, Hac_option)
% EST_2SLS_PILATE Calcola lo stimatore 2SLS sul sottocampione ottimo identificato
% tramite la procedura Bai-Perron (Phase 1).
%
% INPUTS:
%   Y, D, Z, X : Vettori/Matrici dei dati originali
%   pi0        : Frazione del campione
%   m          : Numero di break strutturali (es. 5)
%   Hac_option : (Opzionale) Booleano. Se true, applica Newey-West HAC.
%                Se false o omesso, applica White HC1.

    % --- BACKWARD COMPATIBILITY BLOCK ---
    % Se invocato dalle simulazioni di Sezione 5 (6 argomenti), forza HC1
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
    
    % Separazione algebrica rigorosa: Regressori Proiettati vs Regressori Osservati
    Regressors_proj = [D_hat_c, X_c];
    Regressors_obs  = [D_c, X_c];
    
    % Passaggio del parametro HAC alla subroutine di calcolo della covarianza
    [beta_hat_struct, SE_robust_struct, ~, ~, t_stat_struct] = OLS_Robust_2SLS(Y_c, Regressors_proj, Regressors_obs, T_eff, Hac_option);
    
    beta_val = beta_hat_struct(1);
    SE_val = SE_robust_struct(1);
    t_stat_val = t_stat_struct(1);
end

function [beta_hat, SE_robust, CI, resid, t_stat] = OLS_Robust_2SLS(Y, X_proj, X_obs, T_eff, Hac_option)
% Calcola lo stimatore 2SLS e l'inferenza asintotica (HAC o White-robusta).
% Impiega X_proj per la stima puntuale e X_obs per i residui strutturali.

    [~, nvar] = size(X_proj);

    if size(Y, 2) > 1, Y = Y'; end

    % 1. Stima Puntuale (Invariante rispetto all'HAC)
    beta_hat = (X_proj' * X_proj) \ (X_proj' * Y);

    % 2. Residui Strutturali Corretti (Basati sulle endogene osservate)
    resid = Y - X_obs * beta_hat;

    % 3. Varianza Asintotica
    XX_inv = (X_proj' * X_proj) \ eye(nvar); 
    
    if Hac_option
        Lag = floor(T_eff^(1/3));
        % NeweyWest calcola la somma dei cross-prodotti e divide gia' per (T_eff - nvar).
        % Quindi J_hac = Q_sum / (T_eff - nvar).
        J_hac = NeweyWest(resid, X_proj, Lag, 0);
        
        % La formula sandwich HC1 è: (T_eff / (T_eff - nvar)) * XX_inv * Q_sum * XX_inv
        % Sostituendo Q_sum = J_hac * (T_eff - nvar), i gradi di liberta' si semplificano.
        % Il fattore di scala risultante che moltiplica la struttura diventa esattamente T_eff.
        VCE_robust = T_eff * (XX_inv * J_hac * XX_inv);
    else
        % --- LEGACY BLOCK: Riproduzione esatta al bit per Sezione 5 (HC1 White) ---
        X_res = X_proj .* resid; 
        meat = X_res' * X_res;
        scale = T_eff / (T_eff - nvar);
        VCE_robust = scale * (XX_inv * meat * XX_inv);
    end

    % 4. Inferenza
    SE_robust = sqrt(diag(VCE_robust));
    t_stat = beta_hat ./ SE_robust;
    CI = [beta_hat - 1.96 * SE_robust, beta_hat + 1.96 * SE_robust];
end