
% ========================================================================
% HELPER FUNCTION: computeMMtests (Full Implementation)
%
% CORRECTED: This version correctly "partials out" the effect of the
% exogenous variable X_raw_s_iter from the outcome (Y), endogenous (D),
% and instrument (Z) variables before constructing the moment conditions.
% ========================================================================
function stats_MM = computeMMtests2(Y_s_iter, D_s_iter, Z_raw_s_iter, X_raw_s_iter, T, beta_H0, params_MM)
    % Unpack parameters
    hac_method_MM = params_MM.hac_method_MM;
    trim_prop_cv_mm = params_MM.trim_prop_cv_mm;
    q_instr_MM = params_MM.q_instr_MM;
    Mx_s_local_full_T = params_MM.Mx_s_local_full_T; % Projection matrix for X
    Lag_NW_s_local = params_MM.Lag_NW_s_local;
    test_names_MM_fields = params_MM.test_names_MM_fields;

    stats_MM = struct();
    for mm_test_idx = 1:length(test_names_MM_fields)
        stats_MM.(test_names_MM_fields{mm_test_idx}) = NaN;
    end
    
    % --- CORRECTED SECTION START ---
    % Purge the influence of X from all model variables. This is the
    % standard procedure for handling included exogenous regressors.
    
    % Purge instruments
    Z_instr_MM_s_local = Mx_s_local_full_T*Z_raw_s_iter;

    % Purge outcome and endogenous variables
    Y_purged = Mx_s_local_full_T*Y_s_iter;
    D_purged = Mx_s_local_full_T*D_s_iter;
    
    % Calculate the structural residual using the PURGED variables
    u_structural_s_local = Y_purged - D_purged*beta_H0;
    
    % The Jacobian also needs to be calculated from purged data
    q_t_iter_MM = -Z_instr_MM_s_local .* D_purged;
    % --- CORRECTED SECTION END ---

    valid_MM_instr=rank(Z_instr_MM_s_local'*Z_instr_MM_s_local)==q_instr_MM && ~any(all(abs(Z_instr_MM_s_local)<1e-9,1));
    if ~valid_MM_instr, return; end

    V_ff_MM_hat_s_local=[];
    if strcmp(hac_method_MM,'NW'),V_ff_MM_hat_s_local=NeweyWest(u_structural_s_local,Z_instr_MM_s_local,Lag_NW_s_local,0);
    elseif strcmp(hac_method_MM,'SDHAC'),V_ff_MM_hat_s_local=HAC_Sigma_hat(u_structural_s_local,Z_instr_MM_s_local,1,q_instr_MM,T);
    end
    
    valid_Vff_MM=~isempty(V_ff_MM_hat_s_local)&&rcond(V_ff_MM_hat_s_local)>1e-14&&~(isscalar(V_ff_MM_hat_s_local)&&V_ff_MM_hat_s_local<=1e-12);
    if ~valid_Vff_MM, return; end
    
    V_ff_MM_hat_inv_s_local = inv(V_ff_MM_hat_s_local);
    moments_ft_MM_local = Z_instr_MM_s_local .* u_structural_s_local;
    F_T_MM_s_local = sum(moments_ft_MM_local,1)';
    stats_MM.GMM_S_MM = (1/T) * F_T_MM_s_local' * V_ff_MM_hat_inv_s_local * F_T_MM_s_local;

    tl_mm_s_local=ceil(trim_prop_cv_mm*T); tu_mm_s_local=floor((1-trim_prop_cv_mm)*T);
    num_splits_s_local = 0; if tu_mm_s_local >= tl_mm_s_local, num_splits_s_local=tu_mm_s_local-tl_mm_s_local+1; end
    seq_S_theta_s_MM_local = NaN(1,num_splits_s_local);
    if num_splits_s_local > 0
        s_idx_mm_s_local=0;
        for s_break_s_local=tl_mm_s_local:tu_mm_s_local
            s_idx_mm_s_local=s_idx_mm_s_local+1;
            F_s1_T=sum(moments_ft_MM_local(1:s_break_s_local,:),1)'; F_s2_T=sum(moments_ft_MM_local(s_break_s_local+1:T,:),1)';
            T1_s=s_break_s_local; T2_s=T-s_break_s_local;
            if T1_s>0&&T2_s>0, term1 = NaN; term2 = NaN;
                if T1_s > 1e-6, term1=(1/T1_s)*F_s1_T'*V_ff_MM_hat_inv_s_local*F_s1_T; end
                if T2_s > 1e-6, term2=(1/T2_s)*F_s2_T'*V_ff_MM_hat_inv_s_local*F_s2_T; end
                if ~isnan(term1) && ~isnan(term2), seq_S_theta_s_MM_local(s_idx_mm_s_local)=term1+term2; else, seq_S_theta_s_MM_local(s_idx_mm_s_local) = NaN; end
            else, seq_S_theta_s_MM_local(s_idx_mm_s_local) = NaN; end
        end
    end

    V_qq_hat_feasible=NaN; V_fq_hat_feasible=NaN; V_qq_f_hat_feasible=NaN; V_qq_f_hat_inv_feasible=NaN;
    if strcmp(hac_method_MM,'NW')
        V_qq_hat_feasible = NeweyWest(ones(T,1), q_t_iter_MM, Lag_NW_s_local, 0); 
        V_fq_hat_feasible = NeweyWest_N1N2(moments_ft_MM_local, q_t_iter_MM, Lag_NW_s_local, 0);
    elseif strcmp(hac_method_MM,'SDHAC') % Assuming p_param_dim_MM = 1
        stacked_fq_series=[moments_ft_MM_local, q_t_iter_MM]; 
        V_full_fq_hac=HAC_Sigma_hat(ones(T,1),stacked_fq_series,1,size(stacked_fq_series,2),T); 
        if ~any(isnan(V_full_fq_hac(:)))&&rcond(V_full_fq_hac)>1e-14, V_fq_hat_feasible=V_full_fq_hac(1:q_instr_MM,q_instr_MM+1:end); V_qq_hat_feasible=V_full_fq_hac(q_instr_MM+1:end,q_instr_MM+1:end); end
    end
    if ~any(isnan(V_qq_hat_feasible(:))) && ~any(isnan(V_fq_hat_feasible(:))) && ~any(isnan(V_ff_MM_hat_s_local(:)))
        V_qf_hat_feasible=V_fq_hat_feasible'; 
        if rcond(V_ff_MM_hat_s_local) > 1e-14
            V_qq_f_hat_feasible = V_qq_hat_feasible - V_qf_hat_feasible * V_ff_MM_hat_inv_s_local * V_fq_hat_feasible;
            if ~any(isnan(V_qq_f_hat_feasible(:))) && all(diag(V_qq_f_hat_feasible) > 1e-9) && rcond(V_qq_f_hat_feasible)>1e-14, V_qq_f_hat_inv_feasible = inv(V_qq_f_hat_feasible); end
        end
    end
    if ~any(isnan(V_qq_f_hat_inv_feasible(:))) && num_splits_s_local > 0
        seq_break_criterion_feasible = NaN(1, num_splits_s_local); s_idx_break_sel = 0;
        for s_break_val = tl_mm_s_local:tu_mm_s_local
            s_idx_break_sel=s_idx_break_sel+1; T1_s_val=s_break_val; T2_s_val=T-s_break_val; if T1_s_val==0||T2_s_val==0, continue; end
            F_T1_sel=sum(moments_ft_MM_local(1:T1_s_val,:),1)'; F_T2_sel=sum(moments_ft_MM_local(T1_s_val+1:T,:),1)';
            Q_T1_sel=sum(q_t_iter_MM(1:T1_s_val,:),1)'; Q_T2_sel=sum(q_t_iter_MM(T1_s_val+1:T,:),1)';
            D_T1_feas=Q_T1_sel-V_qf_hat_feasible*V_ff_MM_hat_inv_s_local*F_T1_sel; D_T2_feas=Q_T2_sel-V_qf_hat_feasible*V_ff_MM_hat_inv_s_local*F_T2_sel;
            crit_val1_sel_scaled = NaN; crit_val2_sel_scaled = NaN;
            if T1_s_val > 1e-6, crit_val1_sel_scaled = (T/T1_s_val)*D_T1_feas'*V_qq_f_hat_inv_feasible*D_T1_feas; end 
            if T2_s_val > 1e-6, crit_val2_sel_scaled = (T/T2_s_val)*D_T2_feas'*V_qq_f_hat_inv_feasible*D_T2_feas; end 
            if ~isnan(crit_val1_sel_scaled) && ~isnan(crit_val2_sel_scaled), seq_break_criterion_feasible(s_idx_break_sel)=crit_val1_sel_scaled+crit_val2_sel_scaled; else, seq_break_criterion_feasible(s_idx_break_sel)=NaN; end
        end
        valid_crit_feasible=seq_break_criterion_feasible(~isnan(seq_break_criterion_feasible));
        if ~isempty(valid_crit_feasible), [~,max_crit_idx_s]=max(valid_crit_feasible); original_indices_feasible=find(~isnan(seq_break_criterion_feasible));
            if ~isempty(original_indices_feasible)&&max_crit_idx_s<=length(original_indices_feasible), rt_mm_idx_feasible=original_indices_feasible(max_crit_idx_s);
                if rt_mm_idx_feasible<=length(seq_S_theta_s_MM_local)&&~isnan(seq_S_theta_s_MM_local(rt_mm_idx_feasible)),stats_MM.split_S_MM_Feasible=seq_S_theta_s_MM_local(rt_mm_idx_feasible);end
            end
        end
    end

    if ~isempty(seq_S_theta_s_MM_local)
        valid_S_splits=seq_S_theta_s_MM_local(~isnan(seq_S_theta_s_MM_local));
        if ~isempty(valid_S_splits)
            stats_MM.ave_S_MM=mean(valid_S_splits); 
            S_T_null_for_exp_s_local=stats_MM.GMM_S_MM; X_vals_for_exp_s_local=valid_S_splits-S_T_null_for_exp_s_local; 
            if ~isempty(X_vals_for_exp_s_local) && ~any(isnan(X_vals_for_exp_s_local)) && ~isnan(S_T_null_for_exp_s_local)
                half_X_vals = 0.5 * X_vals_for_exp_s_local; max_half_X = max(half_X_vals); log_mean_exp_half_Stilde_val_s_local = NaN;
                if isinf(max_half_X), log_mean_exp_half_Stilde_val_s_local = Inf; 
                else, mean_exp_half_Stilde = mean(exp(half_X_vals - max_half_X));
                    if mean_exp_half_Stilde > 0, log_mean_exp_half_Stilde_val_s_local = max_half_X + log(mean_exp_half_Stilde); else, log_mean_exp_half_Stilde_val_s_local = -Inf; end
                end
                if ~isinf(log_mean_exp_half_Stilde_val_s_local) && ~isnan(log_mean_exp_half_Stilde_val_s_local), stats_MM.exp_S_MM = 2 * log_mean_exp_half_Stilde_val_s_local + S_T_null_for_exp_s_local; end
            end
        end
    end
    
    if ~isnan(stats_MM.GMM_S_MM) && rcond(V_ff_MM_hat_s_local) > 1e-14
        c_10_qll=10;term_c_qll=c_10_qll/(1+c_10_qll);rtilde_qll=1-c_10_qll/T;
        V_ff_inv_sqrt_qLL=V_ff_MM_hat_s_local^(-0.5);v_t_std_qLL=moments_ft_MM_local*V_ff_inv_sqrt_qLL; qLL_Stilde_sum_NW=0;
        for dim_idx_qll=1:q_instr_MM
            v_t_curr_dim=v_t_std_qLL(:,dim_idx_qll);v_t_curr_dim_demeaned=v_t_curr_dim-mean(v_t_curr_dim); TSSR_N_curr_dim=sum(v_t_curr_dim_demeaned.^2);delta_v_curr_dim=[v_t_curr_dim(1);diff(v_t_curr_dim)]; 
            w_curr_dim=filter(1,[1,-rtilde_qll],delta_v_curr_dim); r_T_vec_proj=rtilde_qll.^((1:T)'); 
            if isempty(r_T_vec_proj)||(r_T_vec_proj'*r_T_vec_proj)<1e-12||rcond(r_T_vec_proj'*r_T_vec_proj)<1e-14,M_vr_dim_proj=eye(T);else,M_vr_dim_proj=eye(T)-r_T_vec_proj*((r_T_vec_proj'*r_T_vec_proj)\r_T_vec_proj');end
            residuals_G_curr_dim=M_vr_dim_proj*w_curr_dim; TSSR_G_curr_dim=sum(residuals_G_curr_dim.^2);
            qLL_Stilde_sum_NW=qLL_Stilde_sum_NW+(TSSR_N_curr_dim-rtilde_qll*TSSR_G_curr_dim);
        end
        stats_MM.qLL_S_MM=term_c_qll*stats_MM.GMM_S_MM+qLL_Stilde_sum_NW;
    end
end % End of computeMMtests
