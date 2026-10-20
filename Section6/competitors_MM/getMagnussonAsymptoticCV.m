function cv = getMagnussonAsymptoticCV(testName, k_moment_conds, alpha_level, N_cv_sims, T_cv_approx, trim_prop, p_zeta)
% getMagnussonAsymptoticCV Simulates asymptotic critical values for qLL-S, exp-S, ave-S.
% This function is based on simulating the limiting distributions from M&M Theorems 1, 2 & 6.
%
% For exp-S, this version implements the form:
%   log(mean(exp(stability_component))) + full_sample_component
% to test against published tables.
%
% Inputs:
%   testName: String, one of 'qLL-S', 'exp-S', 'ave-S'
%   k_moment_conds: Integer, number of moment conditions (k in M&M tables)
%   alpha_level: Scalar, significance level (e.g., 0.05)
%   N_cv_sims: Integer, number of Monte Carlo simulations for CV generation
%   T_cv_approx: Integer, number of points to approximate Brownian motion
%   trim_prop: Scalar, trimming proportion for exp-S and ave-S
%   p_zeta: (Optional) Integer, number of strongly identified parameters concentrated out.
%           Defaults to 0 if not provided.
%
% Output:
%   cv: Scalar, the critical value

    % Set default for p_zeta if not provided
    if nargin < 7
        p_zeta = 0; % Default: no parameters concentrated out
    end

    ds_cv = 1/T_cv_approx;
    % --- Generate Basic Stochastic Components ---
    % 1. Standard k-dimensional Brownian Motion Paths (Wiener Process)
    %    The dimensionality of the Brownian motion is k_moment_conds (k)
    W_all_cv_sims = zeros(T_cv_approx + 1, k_moment_conds, N_cv_sims);
    for n_bm = 1:N_cv_sims
        W_all_cv_sims(:,:,n_bm) = [zeros(1, k_moment_conds); ...
                                 cumsum(randn(T_cv_approx, k_moment_conds) * sqrt(ds_cv), 1)];
    end

    % 2. Chi-squared draws for the full-sample component
    %    Degrees of freedom should be k - p_zeta
    df_chi2_cv = k_moment_conds - p_zeta;

    if df_chi2_cv < 0
        error('CV Sim: df_chi2_cv (k - p_zeta) is negative. k_moment_conds must be >= p_zeta.');
    elseif df_chi2_cv == 0
        chi2_cv_draws = zeros(N_cv_sims, 1);
    else
        chi2_cv_draws = chi2rnd(df_chi2_cv, N_cv_sims, 1);
    end

    test_stats_cv = zeros(N_cv_sims, 1);
    % --- Calculate Test Statistics for each Monte Carlo Replication ---
    for n = 1:N_cv_sims
        W_cv = W_all_cv_sims(:,:,n);
        current_chi2_cv_draw = chi2_cv_draws(n);

        if strcmp(testName, 'qLL-S')
            c_qLL_cv = 10; % Stability component parameter from M&M recommendation [cite: 455]
            
            psi_c_sum_dims_local_cv = 0;
            for dim = 1:k_moment_conds % Looping over k dimensions for stability part
                W_dim_cv = W_cv(:, dim);
                J_dim_cv = zeros(T_cv_approx + 1, 1);
                % Discretized Ornstein-Uhlenbeck process: J(s) = W(s) - c*int(exp(-c(s-r))W(r)dr)
                % or dJ(s) = -c*J(s)ds + dW(s) -> J(t) = (1-c*ds)J(t-1) + (W(t)-W(t-1))
                for t_idx = 1:T_cv_approx
                    J_dim_cv(t_idx+1) = (1 - c_qLL_cv * ds_cv) * J_dim_cv(t_idx) + (W_dim_cv(t_idx+1) - W_dim_cv(t_idx));
                end
                J1_val_cv = J_dim_cv(end); % J(1)
                int_J_sq_val_cv = sum(J_dim_cv(1:T_cv_approx).^2) * ds_cv; % int J(s)^2 ds
                s_grid_cv = (0:ds_cv:1-ds_cv)'; % Grid for integration
                integrand_exp_J_cv = exp(-c_qLL_cv * s_grid_cv) .* J_dim_cv(1:T_cv_approx);
                int_exp_J_val_cv = sum(integrand_exp_J_cv) * ds_cv; % int exp(-cs)J(s)ds

                % Components of psi_c from M&M (2014) Theorem 1 [cite: 570]
                term1_cv = c_qLL_cv * J1_val_cv^2;
                term2_cv = c_qLL_cv^2 * int_J_sq_val_cv;
                term3_num_cv = 2 * c_qLL_cv * (exp(-c_qLL_cv) * J1_val_cv + c_qLL_cv * int_exp_J_val_cv)^2;
                term3_den_cv = 1 - exp(-2 * c_qLL_cv);
                if abs(term3_den_cv) < 1e-12, term3_cv = 0; else, term3_cv = term3_num_cv / term3_den_cv; end
                int_J_val_cv = sum(J_dim_cv(1:T_cv_approx))*ds_cv; % int J(s)ds
                term4_cv = (J1_val_cv + c_qLL_cv * int_J_val_cv)^2;

                psi_c_component_cv = term1_cv + term2_cv + term3_cv - term4_cv;
                psi_c_sum_dims_local_cv = psi_c_sum_dims_local_cv + psi_c_component_cv;
            end
            qLL_Stilde_cv = psi_c_sum_dims_local_cv; % This is sum_i psi_c_i = Psi_c

            c_bar_qLL_cv = 10; % Weight for the full-sample S-component (chi2 part) [cite: 455]
            test_stats_cv(n) = qLL_Stilde_cv + (c_bar_qLL_cv / (1 + c_bar_qLL_cv)) * current_chi2_cv_draw;

        elseif strcmp(testName, 'exp-S') || strcmp(testName, 'ave-S')
            start_idx_cv = floor(trim_prop * T_cv_approx);
            if start_idx_cv == 0, start_idx_cv = 1; end % Ensure index is at least 1
            end_idx_cv = floor((1 - trim_prop) * T_cv_approx);
            if end_idx_cv > T_cv_approx, end_idx_cv = T_cv_approx; end % Ensure index is not beyond T_cv_approx

            psi_tilde_k_tau_values_local_cv = [];
            if start_idx_cv <= end_idx_cv && start_idx_cv > 0
                num_trimmed_points_cv = end_idx_cv - start_idx_cv + 1;
                psi_tilde_k_tau_values_local_cv = zeros(num_trimmed_points_cv, 1);
                W1_cv = W_cv(end, :); % W(1), k_moment_conds-dimensional row vector

                idx_counter_cv = 0;
                for t_idx = start_idx_cv:end_idx_cv % t_idx corresponds to approximating sT in F_sT
                    idx_counter_cv = idx_counter_cv + 1;
                    tau_val_cv = (t_idx) * ds_cv; 
                    if tau_val_cv > 1-1e-9, tau_val_cv = 1-1e-9; end 
                    if tau_val_cv < 1e-9 || (1 - tau_val_cv) < 1e-9
                        psi_tilde_k_tau_values_local_cv(idx_counter_cv) = 0; % Or some other suitable value like NaN if prefer to exclude
                        continue;
                    end
                    W_tau_cv = W_cv(t_idx+1, :); % W(s)
                    W_tilde_tau_cv = W_tau_cv - tau_val_cv * W1_cv; % Brownian Bridge: W(s) - sW(1) [cite: 588]
                    psi_tilde_k_tau_values_local_cv(idx_counter_cv) = (W_tilde_tau_cv * W_tilde_tau_cv') / (tau_val_cv * (1 - tau_val_cv)); % This is S_tilde_s or psi_tilde_k_tau [cite: 588]
                end
            end

            if strcmp(testName, 'exp-S')
                if isempty(psi_tilde_k_tau_values_local_cv) || any(isnan(psi_tilde_k_tau_values_local_cv))
                    test_stats_cv(n) = NaN;
                else
                    % MODIFIED SECTION FOR exp-S (Form 1)
                    % Implements: log(mean(exp(stability_component))) + full_sample_component
                    X_vals = psi_tilde_k_tau_values_local_cv;
                    if isempty(X_vals)
                        log_mean_exp_Stilde_form1 = NaN;
                    else
                        maxX = max(X_vals);
                        if isinf(maxX) % if any X_val was already Inf (unlikely for well-behaved psi_tilde_k_tau)
                             log_mean_exp_Stilde_form1 = Inf;
                        else
                             % log-sum-exp trick for numerical stability of log(mean(exp(X)))
                             log_mean_exp_Stilde_form1 = maxX + log(mean(exp(X_vals - maxX)));
                        end
                    end
                    test_stats_cv(n) = log_mean_exp_Stilde_form1 + current_chi2_cv_draw;
                end
            elseif strcmp(testName, 'ave-S')
                 if isempty(psi_tilde_k_tau_values_local_cv) || any(isnan(psi_tilde_k_tau_values_local_cv))
                    test_stats_cv(n) = NaN;
                 else
                    ave_Stilde_cv = mean(psi_tilde_k_tau_values_local_cv); % This is mean(stability_component)
                    test_stats_cv(n) = ave_Stilde_cv + current_chi2_cv_draw; % Correct for ave-S
                 end
            end
        else
            error('Unknown testName for Magnusson CV generation: %s', testName);
        end
    end

    non_nan_stats = test_stats_cv(~isnan(test_stats_cv));
    if isempty(non_nan_stats)
        cv = NaN; 
        warning('All simulated test statistics were NaN for test: %s, k=%d, p_zeta=%d', testName, k_moment_conds, p_zeta);
    else
        cv = quantile(non_nan_stats, 1 - alpha_level);
    end
end