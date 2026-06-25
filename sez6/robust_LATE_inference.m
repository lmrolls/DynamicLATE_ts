function [stats_AMS, stats_ST, cv_LR] = robust_LATE_inference(Y, D, Z, X, beta_H0, params, user_seed)
% ROBUST_LATE_INFERENCE Computes AMS and ST statistics for Dynamic LATE.
% 'Plug & Play' function for empirical datasets or Monte Carlo simulations.

    % --- 1. INPUT HANDLING AND DEFAULTS ---
    if nargin < 7 || isempty(user_seed), user_seed = randi(1000000); end
    if nargin < 6 || isempty(params)
        params = struct('q', size(Z,2), 'max_breaks', 5, 'eps1', 0.05, 'vector_piL', (0.6:0.1:1), 'step_pi', 0.05);
    end
    
    if ~isfield(params, 'Hac_option'), params.Hac_option = true; end
    if ~isfield(params, 'legacy_math'), params.legacy_math = false; end
    
    if nargin < 5 || isempty(beta_H0), beta_H0 = 0.0; end
    if nargin < 4, X = []; end
    
    T = size(Y, 1);
    q = size(Z, 2); 
    p = size(X, 2);
    a_vec = [beta_H0; 1]; 
    b_vec = [1; -beta_H0];

    stats_AMS = struct('AR', NaN, 'LM', NaN, 'LR', NaN, 'M2', NaN);
    stats_ST  = struct('AR', NaN, 'LM', NaN, 'LR', NaN, 'M2', NaN, 'best_pi', NaN, 'best_m', NaN, 'dates', []);
    cv_LR     = struct('cv_AMS', NaN, 'cv_ST', NaN, 'pval_AMS', NaN, 'pval_ST', NaN);

    % --- 2. ORTHOGONAL PROJECTION (M_X) ---
    Mx = eye(T);
    if p > 0 && rank(X'*X) == p
        % [LEGACY GUARANTEE]: Even under legacy_math, the global Mx used mldivide (\) in the original script.
        Mx = eye(T) - X * ((X'*X) \ X'); 
    end
    Zbar_AMS = Mx * Z;

    % --- 3. COMPUTATION OF AMS STATISTICS ---
    res_AMS = computeStats_Internal(Y, D, Zbar_AMS, X, T, q, p, a_vec, b_vec, false, params.Hac_option, params.legacy_math);
    stats_AMS.AR = res_AMS.M1 / q;
    stats_AMS.LM = (res_AMS.M12^2) / res_AMS.M2;
    stats_AMS.LR = 0.5 * (res_AMS.M1 - res_AMS.M2 + sqrt(max(0, (res_AMS.M1 - res_AMS.M2)^2 + 4 * res_AMS.M12^2)));
    stats_AMS.M2 = res_AMS.M2;
    
    % RNG Branching logic to preserve 12/06 paper reproducibility
    if params.legacy_math
        [cv_LR.cv_AMS, cv_LR.pval_AMS] = LR_cv_legacy_exact(stats_AMS.M2, q, 0.05, user_seed + 1, stats_AMS.LR);
    else
        [cv_LR.cv_AMS, cv_LR.pval_AMS] = LR_cv_simulate_internal(stats_AMS.M2, q, 0.05, user_seed + 1, stats_AMS.LR);
    end

    % --- 4. COMPUTATION OF ST STATISTICS ---
    if params.max_breaks > 0 && exist('Procedure_piL_Full_Sample', 'file') == 2
        params.h_proc = round(params.eps1 * T);
        [~, datevec_ST, ~] = Procedure_piL_Full_Sample(T, D, Z, q, params.max_breaks, params.h_proc, params.eps1, 1, 0, 1, 1, 1, 0, X, p, 0.0001, 20, 0, 0, 1, 1);
        
        Sup_M2_piL = NaN(length(params.vector_piL), 1); 
        Sup_m_piL = NaN(length(params.vector_piL), 1); 
        Sup_pi_piL = NaN(length(params.vector_piL), 1);
        
        for piL_idx = 1:length(params.vector_piL)
            num_pis = params.vector_piL(piL_idx):params.step_pi:1; 
            M2_pi_m = NaN(length(num_pis), params.max_breaks);
            for m_idx = 1:params.max_breaks
                for pi_idx = 1:length(num_pis)
                    Index_S = zeros(T, 1);
                    if m_idx == 0
                        Index_S = ones(T, 1);
                    elseif m_idx <= size(datevec_ST, 2) && any(datevec_ST(:, m_idx) ~= 0)
                        S_seg = S_Dates_m_M2_Core(Y, D, Z, X, Mx, beta_H0, datevec_ST, T, num_pis(pi_idx), m_idx, q, p, a_vec, b_vec, params.Hac_option, params.legacy_math);
                        if ~isempty(S_seg)
                            for r = 1:size(S_seg, 1)
                                if S_seg(r,1) >= 1 && S_seg(r,2) <= T, Index_S(S_seg(r,1):S_seg(r,2)) = 1; end
                            end
                        end
                    else
                        continue; 
                    end
                    
                    Z_ST_tmp = Z; Z_ST_tmp(Index_S == 0, :) = 0;
                    if sum(Index_S) > max(q, p)
                        res_tmp = computeStats_Internal(Y, D, Mx * Z_ST_tmp, X, T, q, p, a_vec, b_vec, false, params.Hac_option, params.legacy_math);
                        M2_pi_m(pi_idx, m_idx) = res_tmp.M2;
                    end
                end
            end
            [Val_m, ind_m] = max(M2_pi_m, [], 1); [Val_m_pi, ind_m_pi] = max(Val_m);
            if ~isempty(Val_m_pi) && ~isnan(Val_m_pi)
                Sup_M2_piL(piL_idx) = Val_m_pi; 
                Sup_m_piL(piL_idx) = ind_m_pi; 
                Sup_pi_piL(piL_idx) = num_pis(ind_m(ind_m_pi)); 
            end
        end
        
        [~, idx_Sup] = max(Sup_M2_piL);
        if ~isempty(idx_Sup) && ~isnan(Sup_m_piL(idx_Sup))
            S_seg_final = S_Dates_m_M2_Core(Y, D, Z, X, Mx, beta_H0, datevec_ST, T, Sup_pi_piL(idx_Sup), Sup_m_piL(idx_Sup), q, p, a_vec, b_vec, params.Hac_option, params.legacy_math);
            Index_S_final = zeros(T, 1);
            if ~isempty(S_seg_final)
                for r = 1:size(S_seg_final, 1)
                    if S_seg_final(r,1) >= 1 && S_seg_final(r,2) <= T, Index_S_final(S_seg_final(r,1):S_seg_final(r,2)) = 1; end
                end
            end
            Z_ST_final = Z; Z_ST_final(Index_S_final == 0, :) = 0;
            if sum(Index_S_final) > max(q, p)
                res_ST_final = computeStats_Internal(Y, D, Mx * Z_ST_final, X, T, q, p, a_vec, b_vec, true, params.Hac_option, params.legacy_math);
                stats_ST.AR = res_ST_final.M1 / q;
                stats_ST.LM = (res_ST_final.M12^2) / res_ST_final.M2;
                stats_ST.LR = 0.5 * (res_ST_final.M1 - res_ST_final.M2 + sqrt(max(0, (res_ST_final.M1 - res_ST_final.M2)^2 + 4 * res_ST_final.M12^2)));
                stats_ST.M2 = res_ST_final.M2;
                stats_ST.best_pi = Sup_pi_piL(idx_Sup);
                stats_ST.best_m = Sup_m_piL(idx_Sup);
                stats_ST.dates = S_seg_final;
                
                if params.legacy_math
                    [cv_LR.cv_ST, cv_LR.pval_ST] = LR_cv_legacy_exact(stats_ST.M2, q, 0.05, user_seed + 2, stats_ST.LR);
                else
                    [cv_LR.cv_ST, cv_LR.pval_ST] = LR_cv_simulate_internal(stats_ST.M2, q, 0.05, user_seed + 2, stats_ST.LR);
                end
            end
        end
    end
end

function [cv, pval] = LR_cv_legacy_exact(M2_ST, q, alpha, seed, LR_stat)
    % DE-VECTORIZED: Replicates the exact RNG extraction order of the 12/06 legacy script
    rng(seed, "twister");
    n_reps = 10000;
    LR_ST = zeros(n_reps, 1);
    sqrt_M2 = sqrt(M2_ST);
    for i_reps = 1:n_reps
        M1_ST = chi2rnd(q);
        M12_ST = sqrt_M2 * randn;
        LR_ST(i_reps, 1) = 0.5 * (M1_ST - M2_ST + sqrt((M1_ST - M2_ST)^2 + 4 * M12_ST^2));
    end
    cv = quantile(LR_ST, 1 - alpha);
    pval = mean(LR_ST >= LR_stat); % Extract exact p-value from ECDF
end

function [cv, pval] = LR_cv_simulate_internal(M2_val, q, alpha, s_seed, LR_stat)
    % MODERN VECTORIZED EXTRACTION
    if isnan(M2_val), cv = NaN; pval = NaN; return; end
    rng(s_seed, "twister");
    num_reps = 10000;
    N1_sim = randn(num_reps, q);
    N2_sim = randn(num_reps, q) + sqrt(M2_val);
    LR_sim = zeros(num_reps, 1);
    for i = 1:num_reps
        m1 = N1_sim(i,:) * N1_sim(i,:)';
        m2 = N2_sim(i,:) * N2_sim(i,:)';
        m12 = N1_sim(i,:) * N2_sim(i,:)';
        LR_sim(i) = 0.5 * (m1 - m2 + sqrt((m1 - m2)^2 + 4 * m12^2));
    end
    cv = quantile(LR_sim, 1 - alpha);
    pval = mean(LR_sim >= LR_stat); % Extract exact p-value from ECDF
end

function stats = computeStats_Internal(Y, D, Zbar, X, T_eff, q, p, a, b, is_ST, Hac_option, legacy_math)
    y = [Y, D]; 
    if isempty(X)
        if legacy_math
            v_hat = y - Zbar * inv(Zbar' * Zbar) * Zbar' * y;
        else
            v_hat = y - Zbar * ((Zbar' * Zbar) \ Zbar') * y;
        end
    else
        if legacy_math
            v_hat = y - Zbar * inv(Zbar' * Zbar) * Zbar' * y - X * inv(X' * X) * X' * y; 
        else
            v_hat = y - Zbar * ((Zbar' * Zbar) \ Zbar') * y - X * ((X' * X) \ X') * y; 
        end
    end
    Sigma_v_hat = (v_hat' * v_hat) / (T_eff - q - p); 
    
    if Hac_option, Lag = floor(T_eff^(1/3)); else, Lag = 0; end
    
    Sigma_N1_hat_main = NeweyWest(v_hat * b, Zbar, Lag, 0);
    Sigma_N2star_hat  = NeweyWest(v_hat * inv(Sigma_v_hat) * a, Zbar, Lag, 0);
    Sigma_N1N2_hat    = NeweyWest_N1N2((v_hat * b .* Zbar), (v_hat * inv(Sigma_v_hat) * a .* Zbar), Lag, 0);
    Sigma_N2_hat      = Sigma_N2star_hat - Sigma_N1N2_hat * inv(Sigma_N1_hat_main) * Sigma_N1N2_hat';
    
    if is_ST
        Sigma_N1_hat_lag0 = NeweyWest(v_hat * b, Zbar, 0, 0);
        N1_val = min(Sigma_N1_hat_main^(-1/2), Sigma_N1_hat_lag0^(-1/2)) * T_eff^(-1/2) * Zbar' * y * b;
    else
        N1_val = Sigma_N1_hat_main^(-1/2) * T_eff^(-1/2) * Zbar' * y * b; 
    end
    N2_val = Sigma_N2_hat^(-1/2) * (T_eff^(-1/2) * Zbar' * y * inv(Sigma_v_hat) * a - Sigma_N1N2_hat * (Sigma_N1_hat_main^(-1/2)) * N1_val);
    
    stats.M1 = N1_val' * N1_val; 
    stats.M12 = N1_val' * N2_val; 
    stats.M2 = N2_val' * N2_val;
end

function [S_m] = S_Dates_m_M2_Core(Y, D, Z, X, Mx, beta_H0, datevec, T, pi, m, q, p, a_vec, b_vec, Hac_option, legacy_math)
    Tpi = floor(pi * T); vec_M2_ST = zeros(T, 1);
    if m == 1, Break_Date = [1, datevec(1,1); datevec(1,1)+1, T];
    elseif m == 2, Break_Date = [1, datevec(1,2); datevec(1,2)+1, datevec(2,2); datevec(2,2)+1, T];
    elseif m == 3, Break_Date = [1, datevec(1,3); datevec(1,3)+1, datevec(2,3); datevec(2,3)+1, datevec(3,3); datevec(3,3)+1, T];
    elseif m == 4, Break_Date = [1, datevec(1,4); datevec(1,4)+1, datevec(2,4); datevec(2,4)+1, datevec(3,4); datevec(3,4)+1, datevec(4,4); datevec(4,4)+1, T];
    else, Break_Date = [1, datevec(1,5); datevec(1,5)+1, datevec(2,5); datevec(2,5)+1, datevec(3,5); datevec(3,5)+1, datevec(4,5); datevec(4,5)+1, datevec(5,5); datevec(5,5)+1, T]; end
    Table_t_stat_dates = zeros(m+1, 3);
    for i = 1:m+1
        CT_Z_S = Z(Break_Date(i,1):Break_Date(i,2), :); T_S = size(CT_Z_S, 1);
        if isempty(X)
            Mx_S = eye(T_S); X_S = []; 
        else
            X_S = X(Break_Date(i,1):Break_Date(i,2), :); 
            if legacy_math, Mx_S = eye(T_S) - X_S * inv(X_S' * X_S) * X_S'; 
            else, Mx_S = eye(T_S) - X_S * ((X_S' * X_S) \ X_S'); end
        end
        stats = computeStats_Internal(Y(Break_Date(i,1):Break_Date(i,2), 1), D(Break_Date(i,1):Break_Date(i,2), 1), Mx_S * CT_Z_S, X_S, T_S, q, p, a_vec, b_vec, false, Hac_option, legacy_math);
        Table_t_stat_dates(i,1:3) = [stats.M2, Break_Date(i,1), Break_Date(i,2)];
    end
    Table_t_stat_dates = sortrows(Table_t_stat_dates, 1, 'descend');
    A = cumsum(Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2) + 1, 1); Table_t_stat_dates = [Table_t_stat_dates A];
    AA = A >= Tpi; AA = AA & cumsum(AA, 1) < 2;
    [ind_larger, ~] = find(AA); ind_smaller = ind_larger - 1;
    if ind_smaller < 1, Missing = Tpi; Start_date = Table_t_stat_dates(ind_larger, 2) - 1; End_date = Table_t_stat_dates(ind_larger, 3);
    else, Missing = Tpi - A(ind_smaller, 1); Start_date = Table_t_stat_dates(ind_larger, 2) - 1; End_date = Table_t_stat_dates(ind_larger, 3); end
    if (End_date - Start_date) == Missing
        if m == 1, pi_Dates = [1, T];
        elseif ind_smaller < 1, pi_Dates = [max(1, Start_date), End_date-1];
        else, pi_Dates = sortrows([Table_t_stat_dates(1:ind_smaller, 2:3); max(1, Start_date)+1, End_date], 1, 'ascend'); end
    else
        for i = 1:(End_date - Start_date)-Missing
            CT_Z_S = Z(Start_date+i:Start_date+i+Missing, :); T_S = size(CT_Z_S, 1);
            if isempty(X)
                Mx_S = eye(T_S); X_S = []; 
            else
                X_S = X(Start_date+i:Start_date+i+Missing, :); 
                if legacy_math, Mx_S = eye(T_S) - X_S * inv(X_S' * X_S) * X_S'; 
                else, Mx_S = eye(T_S) - X_S * ((X_S' * X_S) \ X_S'); end
            end
            stats = computeStats_Internal(Y(Start_date+i:Start_date+i+Missing, 1), D(Start_date+i:Start_date+i+Missing, 1), Mx_S * CT_Z_S, X_S, T_S, q, p, a_vec, b_vec, false, Hac_option, legacy_math);
            vec_M2_ST(i,1) = stats.M2;
        end
        [~, ind] = max(abs(vec_M2_ST));
        pi_Dates = sortrows([Table_t_stat_dates(1:ind_smaller, 2:3); Start_date+ind, Start_date+ind+Missing-1], 1, 'ascend');
    end
    S_m = pi_Dates;
end