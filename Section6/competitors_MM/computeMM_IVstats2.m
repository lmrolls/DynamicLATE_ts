% ========================================================================
% HELPER FUNCTION: computeMM_IVstats2 (Fully Modified & Corrected)
% ========================================================================
% This function computes the "split-sample" family of tests from Magnusson
% & Mavroeidis (2014), specifically for the linear IV case.
%
% MODIFIED: This version now correctly handles the presence of exogenous
% regressors (X) and allows for either Homoskedastic ("HOM") or
% Newey-West ("NW") HAC variance-covariance estimation. It also fixes an
% error with the non-standard 'vec' function.
%
% UPDATED: Now calculates split_KLM, split_JKLM, and split_CLR for the
% homoskedastic ('HOM') case.
% ========================================================================
function stats_MM_IV = computeMM_IVstats2(Y_s_iter, D_s_iter, Z_raw_s_iter, X_raw_s_iter, T, beta_H0, params)
    
    %% --- 1. Unpack Parameters and Initial Setup ---
    p_params = params.p_param_dim;
    k_instr = params.q_instr;
    trim_prop = params.trim_prop_cv_mm;
    hac_method = params.hac_method_MM; % "HOM" or "NW"
    Lag_NW = params.Lag_NW_s_local;

    test_names_MM_IV_fields = {'split_S_IV', 'split_KLM_IV', 'split_JKLM_IV', 'split_CLR_IV', 'rk_stat_IV'};
    stats_MM_IV = struct();
    for test_idx = 1:length(test_names_MM_IV_fields)
        stats_MM_IV.(test_names_MM_IV_fields{test_idx}) = NaN;
    end

    % Unpack the projection matrix for X.
    Mx_s_local_full_T = params.Mx_s_local_full_T;

    %% --- 2. Partial Out Exogenous Regressors (The Correction) ---
    % This is the crucial step. We remove the linear influence of X from 
    % all other model variables by pre-multiplying with the projection matrix Mx.
    Y_purged = Mx_s_local_full_T * Y_s_iter;
    D_purged = Mx_s_local_full_T * D_s_iter;
    Z_purged = Mx_s_local_full_T * Z_raw_s_iter;

    %% --- 3. Main Logic Branching on HAC Method ---
    if strcmp(hac_method, 'HOM')
        % --- 3a. HOMOSKEDASTIC CASE ---
        % This section follows the simplified logic from Section 3.3.1 of the
        % Magnusson & Mavroeidis (2014) paper, assuming homoskedastic errors.
        Y_data = [Y_purged, D_purged];
        A0 = [beta_H0; eye(p_params)];

        try
            first_stage_residuals = D_purged - Z_purged * ((Z_purged' * Z_purged) \ (Z_purged' * D_purged));
            Omega = cov(first_stage_residuals);
            if rcond(Omega) < 1e-12, return; end
            Omega_inv = inv(Omega);
        catch
            return;
        end

        AOA_term = A0' * Omega_inv * A0;
        if rcond(AOA_term) < 1e-12, return; end
        AOA_inv_sqrt = sqrtm(inv(AOA_term));

        t_range = ceil(trim_prop * T) : floor((1 - trim_prop) * T);
        max_criterion_val = -Inf;
        t_b_hat = -1;

        for t_b = t_range
            Z_tau = [Z_purged(1:t_b, :), zeros(t_b, k_instr); ...
                     zeros(T - t_b, k_instr), Z_purged(t_b+1:end, :)];
            ZtZ = Z_tau' * Z_tau;
            if rcond(ZtZ) < 1e-12, continue; end
            ZtZ_inv_sqrt = sqrtm(inv(ZtZ));
            D_bar_tau = ZtZ_inv_sqrt * Z_tau' * Y_data * Omega_inv * A0 * AOA_inv_sqrt;
            criterion_val = trace(D_bar_tau' * D_bar_tau);
            if criterion_val > max_criterion_val
                max_criterion_val = criterion_val;
                t_b_hat = t_b;
            end
        end

        if t_b_hat == -1, return; end
        stats_MM_IV.rk_stat_IV = max_criterion_val;

        b0 = [1; -beta_H0];
        b0Ob0_term = b0' * Omega * b0;
        if b0Ob0_term < 1e-9, return; end
        b0Ob0_inv_sqrt = (b0Ob0_term)^(-0.5);

        Z_tau_hat = [Z_purged(1:t_b_hat, :), zeros(t_b_hat, k_instr); ...
                     zeros(T - t_b_hat, k_instr), Z_purged(t_b_hat+1:end, :)];
        ZtZ_hat = Z_tau_hat' * Z_tau_hat;
        if rcond(ZtZ_hat) < 1e-12, return; end
        ZtZ_hat_inv_sqrt = sqrtm(inv(ZtZ_hat));

        F_bar_hat = ZtZ_hat_inv_sqrt * Z_tau_hat' * Y_data * b0 * b0Ob0_inv_sqrt;
        D_bar_hat = ZtZ_hat_inv_sqrt * Z_tau_hat' * Y_data * Omega_inv * A0 * AOA_inv_sqrt;
        
        stats_MM_IV.split_S_IV = F_bar_hat' * F_bar_hat;

        % --- UPDATED SECTION START: Compute KLM for HOM case ---
        try
            % Project F_bar onto D_bar to get the KLM statistic
            proj_D_bar_hat = D_bar_hat * inv(D_bar_hat' * D_bar_hat) * D_bar_hat';
            stats_MM_IV.split_KLM_IV = F_bar_hat' * proj_D_bar_hat * F_bar_hat;
        catch
            stats_MM_IV.split_KLM_IV = NaN;
        end
        % --- UPDATED SECTION END ---

    elseif strcmp(hac_method, 'NW')
        % --- 3b. NEWEY-WEST (HAC) CASE ---
        % This section follows the more general logic from Section 3.3.2,
        % using HAC estimators for the variance components.
        
        u_structural = Y_purged - D_purged*beta_H0;
        moments_f = Z_purged .* u_structural;
        jacobian_q = -Z_purged .* D_purged;

        V_ff_hat = NeweyWest(ones(T,1), moments_f, Lag_NW, 0);
        V_qq_hat = NeweyWest(ones(T,1), jacobian_q, Lag_NW, 0);
        V_fq_hat = NeweyWest_N1N2(moments_f, jacobian_q, Lag_NW, 0);

        if rcond(V_ff_hat) < 1e-14, return; end
        V_ff_inv = inv(V_ff_hat);
        V_qf_hat = V_fq_hat';
        V_qq_f_hat = V_qq_hat - V_qf_hat * V_ff_inv * V_fq_hat;
        if rcond(V_qq_f_hat) < 1e-14, return; end
        V_qq_f_inv = inv(V_qq_f_hat);

        t_range = ceil(trim_prop * T) : floor((1 - trim_prop) * T);
        max_criterion_val = -Inf;
        t_b_hat = -1;

        for t_b = t_range
            F_T_1 = sum(moments_f(1:t_b, :), 1)';
            F_T_2 = sum(moments_f(t_b+1:T, :), 1)';
            Q_T_1 = sum(jacobian_q(1:t_b, :), 1)';
            Q_T_2 = sum(jacobian_q(t_b+1:T, :), 1)';
            
            % FIXED: Replaced non-standard vec() with MATLAB's colon operator (:)
            vec_D_T_1 = Q_T_1(:) - V_qf_hat * V_ff_inv * F_T_1;
            vec_D_T_2 = Q_T_2(:) - V_qf_hat * V_ff_inv * F_T_2;
            
            criterion = (1/t_b) * vec_D_T_1' * V_qq_f_inv * vec_D_T_1 + (1/(T-t_b)) * vec_D_T_2' * V_qq_f_inv * vec_D_T_2;
            
            if criterion > max_criterion_val
                max_criterion_val = criterion;
                t_b_hat = t_b;
            end
        end

        if t_b_hat == -1, return; end
        stats_MM_IV.rk_stat_IV = max_criterion_val;
        
        % Recompute components at the estimated break date
        F1 = sum(moments_f(1:t_b_hat, :), 1)';
        F2 = sum(moments_f(t_b_hat+1:T, :), 1)';
        Q1 = sum(jacobian_q(1:t_b_hat, :), 1)';
        
        % FIXED: Replaced non-standard vec() with MATLAB's colon operator (:)
        vec_D1 = Q1(:) - V_qf_hat * V_ff_inv * F1;
        D1_mat = reshape(vec_D1, k_instr, p_params);
        
        % Compute split-S
        stats_MM_IV.split_S_IV = (1/t_b_hat) * (F1' * V_ff_inv * F1) + (1/(T-t_b_hat)) * (F2' * V_ff_inv * F2);
        
        % Compute split-KLM
        try
            V_ff_inv_sqrt = sqrtm(V_ff_inv);
            F_tilde = V_ff_inv_sqrt * [F1, F2];
            D_tilde = V_ff_inv_sqrt * D1_mat; % Jacobian is the same in both regimes under this setup
            
            proj_D_tilde = D_tilde * inv(D_tilde' * D_tilde) * D_tilde';
            
            KLM1 = F_tilde(:,1)' * proj_D_tilde * F_tilde(:,1);
            KLM2 = F_tilde(:,2)' * proj_D_tilde * F_tilde(:,2);
            stats_MM_IV.split_KLM_IV = (1/t_b_hat)*KLM1 + (1/(T-t_b_hat))*KLM2;
        catch
             stats_MM_IV.split_KLM_IV = NaN;
        end
    else
        error('Invalid HAC method specified. Must be "HOM" or "NW".');
    end

    %% --- 4. Compute Final Statistics (JKLM, CLR) ---
    % This part is common to both HOM and NW cases
    if ~isnan(stats_MM_IV.split_S_IV) && ~isnan(stats_MM_IV.split_KLM_IV)
        stats_MM_IV.split_JKLM_IV = stats_MM_IV.split_S_IV - stats_MM_IV.split_KLM_IV;
    else
        stats_MM_IV.split_JKLM_IV = NaN;
    end
    
    if ~isnan(stats_MM_IV.split_S_IV) && ~isnan(stats_MM_IV.split_JKLM_IV) && ~isnan(stats_MM_IV.rk_stat_IV)
        s_plus_rk = stats_MM_IV.split_S_IV + stats_MM_IV.rk_stat_IV;
        s_minus_rk = stats_MM_IV.split_S_IV - stats_MM_IV.rk_stat_IV;
        sqrt_term = s_plus_rk^2 - 4 * stats_MM_IV.split_JKLM_IV * stats_MM_IV.rk_stat_IV;

        if sqrt_term >= 0
             stats_MM_IV.split_CLR_IV = 0.5 * (s_minus_rk + sqrt(sqrt_term));
        else
             stats_MM_IV.split_CLR_IV = NaN;
        end
    else
        stats_MM_IV.split_CLR_IV = NaN;
    end

end % End of computeMM_IVstats2 function
