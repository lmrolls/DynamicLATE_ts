function cv_split_CLR = split_CLR_cv_seed(rk_stat, p, k, alpha, seed)
% ========================================================================
% HELPER FUNCTION: split_CLR_cv_seed
% ========================================================================
% Computes the data-dependent critical value for the split-sample CLR test
% statistic, based on Theorem 3, part 5 from Magnusson & Mavroeidis (2014).
%
% Inputs:
%   rk_stat - Scalar, the rank statistic (D_bar' * D_bar) from computeMM_IVstats.
%   p       - Scalar, number of endogenous regressors.
%   k       - Scalar, number of instruments in one regime. The total number
%             of split-sample instruments is 2*k.
%   alpha   - Scalar, significance level (e.g., 0.05).
%   seed    - Scalar, for reproducible random number generation.
%
% Output:
%   cv_split_CLR - Scalar, (1-alpha) quantile of the simulated distribution.
% ========================================================================

    % Set a reproducible random seed
    rng(seed, "twister");

    % Number of simulations to generate the distribution
    n_reps = 10000;
    split_CLR_dist = zeros(n_reps, 1);

    % Define degrees of freedom based on Theorem 3
    df_p = p;
    df_2k_minus_p = 2 * k - p;

    % Check for valid degrees of freedom
    if df_2k_minus_p <= 0
        cv_split_CLR = inf;
        return;
    end

    for i_rep = 1:n_reps
        % Simulate the two independent chi-squared random variables
        psi_p = chi2rnd(df_p);
        psi_2k_minus_p = chi2rnd(df_2k_minus_p);
        
        % The distribution of the split-CLR statistic conditional on rk_stat
        term1 = psi_p + psi_2k_minus_p + rk_stat;
        term2 = psi_p + psi_2k_minus_p - rk_stat;
        sqrt_term = term1^2 - 4 * psi_2k_minus_p * rk_stat;
        
        if sqrt_term >= 0
            split_CLR_dist(i_rep) = 0.5 * (term2 + sqrt(sqrt_term));
        else
            % This case should be rare, but handle it to avoid complex numbers
            split_CLR_dist(i_rep) = 0.5 * term2;
        end
    end

    % The critical value is the (1-alpha) quantile of the simulated distribution
    cv_split_CLR = quantile(split_CLR_dist, 1 - alpha);
end
