function Results = runSection7_Table3(alpha_grid, pairs, verbose)
% RUNSECTION7_TABLE3  Table 3 of the paper (Section 7.3): weak identification-robust inference on beta.
%
% Legacy script: pi-LATE\Identification_Robust_Inference\Run_Identi_Robus_Tests_NS18_Heter.m (8/6/2025), run
% once per (Y, D, alpha) by editing lines 31-33 and 80 (README of the package). Here the script is a function:
% data, sample, construction of the instrument, dating on the control sample, selection of S_hat
% (S_Dates_m_M2_Heter.m, legacy_IRI\), statistics and simulated LR critical values are computed exactly as
% in the script; only the loop over the four published (Y, D) pairs and the three significance levels is new.
%
% SYNTAX:
%   Results = runSection7_Table3()                      % published configurations, alpha in {0.10, 0.05, 0.01}
%   Results = runSection7_Table3([0.10 0.05 0.01], {'DNF2','DNY2'})   % one pair {Y_name, D_name}
%
% DESIGN (Section 7, heteroskedasticity-based identification):
%   * Data: ..\m_break\dataForMatlabAltEstTick.csv, observations 401:end (T = 836); FOMC days (FOMCused = 1,
%     T_P) and control days (FOMCused = 0, T_C). Sample re-ordered as [control; treatment].
%   * D_t: 1-day change in 2-Year Treasury yields (DNY2) or 30-minute policy news (path_intra_wide), demeaned
%     within group; Y_t: 2-Year nominal (DNF2) or real (DRF2) forward rates, not demeaned (as in the script).
%   * Instrument Z_t = (T/T_P) D_t on FOMC days and -(T/T_C) D_t on control days (demeaned D), q = 1, no X.
%   * S_hat: Bai-Perron dating of the control sample only (Procedure_piL_Full_Sample on D(1:T_C), Z(1:T_C),
%     m_breaks = 5, eps1 = 0.05, h = round(eps1*T) with the FULL T, robust = 1), selection over
%     pi_L in 0.7:0.1:1 (not 0.6 as in the simulations), pi in pi_L:0.05:1 and m = 1..5 of the sub-sample
%     maximising M2 computed on the full re-ordered sample with the instrument zeroed outside S
%     (S_Dates_m_M2_Heter: the treatment block [T_C+1, T] is always included, M2 = inf).
%   * Statistics: AR = M1/q, LM = M12^2/M2, LR (CLR) with the NW HAC covariances (lag floor(T^(1/3))), robust
%     and homoskedastic versions, on S_hat (ST) and on the full sample (AMS). Rejection: AR, LM vs chi2 CVs;
%     LR vs the simulated conditional CV LR_critical_value(M2, q, alpha) (1000 draws) after rng(11), in the
%     legacy order ST, ST_hom, AMS, AMS_hom; the three alpha levels use the same draws (quantiles).
%   * The published table reports LM_T, CLR_T (full sample) and LM_T(S_hat), CLR_T(S_hat): a tick = rejection.
%   * The script's load('pi.mat') (an unrelated 0/1 vector named pi) has no effect on the computations and is
%     omitted; the control/treatment dates are not needed (T_C, T_P from FOMCused).
% =========================================================================

    if nargin < 1 || isempty(alpha_grid), alpha_grid = [0.10, 0.05, 0.01]; end
    if nargin < 2 || isempty(pairs)
        pairs = {{'DNF2', 'path_intra_wide'}, {'DRF2', 'path_intra_wide'}, {'DNF2', 'DNY2'}, {'DRF2', 'DNY2'}};
    elseif ~iscell(pairs{1})
        pairs = {pairs};
    end
    if nargin < 3 || isempty(verbose), verbose = true; end

    restoredefaultpath;
    addpath('..\m_break');          % Procedure_piL_Full_Sample (Bai-Perron dating), NeweyWest, NeweyWest_N1N2
    addpath('legacy_IRI');          % S_Dates_m_M2_Heter, LR_critical_value (legacy, verbatim)
    warning('off', 'all');

    folder_name = 'Paper_Exhibits_Results_SEC7_TABLE3';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

    % --- Data (legacy lines 18-28) ---
    dta = readtable(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv'));
    FOMC = dta.FOMCused(401:end);
    T   = numel(FOMC);
    T_C = sum(FOMC == 0);
    T_P = sum(FOMC == 1);

    % --- Legacy options of the dating procedure (lines 84-154) ---
    beta_H0 = 0; m_breaks = 5; eps1 = 0.05; h = round(eps1 * T);
    fixb = 0; betaini = 0; maxi = 20; printd = 0; eps = 0.0001;
    robust = 1; prewhit = 0; hetdat = 1; hetvar = 1; hetomega = 1; hetq = 1; doglobal = 1;
    vector_piL = 0.7:0.1:1; step_pi = 0.05;
    X = []; p = 0; q = 1; Mx = eye(T);
    a = [beta_H0; 1]; b = [1; -beta_H0];
    Lag = floor(T^(1/3));

    Results = struct();
    Results.alpha_grid = alpha_grid; Results.T = T; Results.T_C = T_C; Results.T_P = T_P;
    labels = {'AR_ST', 'LM_ST', 'LR_ST', 'AR_ST_hom', 'LM_ST_hom', 'LR_ST_hom', 'AR_AMS', 'LM_AMS', 'LR_AMS', 'AR_AMS_hom', 'LM_AMS_hom', 'LR_AMS_hom'};

    for ip = 1:numel(pairs)
        Y_name = pairs{ip}{1}; D_name = pairs{ip}{2};
        lbl = sprintf('Y_%s_D_%s', Y_name, D_name);
        if verbose, fprintf('\n=== Table 3: Y = %s, D = %s ===\n', Y_name, D_name); end

        % --- Variables and instrument (legacy lines 31-79) ---
        D_raw = dta.(D_name)(401:end); Y_raw = dta.(Y_name)(401:end);
        assert(~any(isnan(D_raw)) && ~any(isnan(Y_raw)), 'NaN in %s or %s over the sample 401:end.', D_name, Y_name);
        Yc = Y_raw(FOMC == 0); Yp = Y_raw(FOMC == 1);
        Y = [Yc; Yp];                                   % not demeaned (legacy: the demeaned version is commented out)
        Dc = D_raw(FOMC == 0); Dp = D_raw(FOMC == 1);
        D = [Dc - mean(Dc, 1); Dp - mean(Dp, 1)];
        Z = [-T/T_C * (Dc - mean(Dc, 1)); T/T_P * (Dp - mean(Dp, 1))];
        Y_orig = Y; D_orig = D; Z_orig = Z;

        % --- Bai-Perron dating on the control sample (legacy line 159) ---
        [~, datevec, ~] = Procedure_piL_Full_Sample(T_C, D(1:T_C, :), Z(1:T_C, :), q, m_breaks, h, eps1, robust, prewhit, ...
                                                    hetomega, hetq, doglobal, fixb, X, p, eps, maxi, betaini, printd, hetdat, hetvar);

        % --- Selection of S_hat (legacy lines 167-257) ---
        n_piL = numel(vector_piL);
        Sup_M2_piL = NaN(n_piL, 1); Sup_m_piL = NaN(n_piL, 1); Sup_pi_piL = NaN(n_piL, 1);
        for piL = 1:n_piL
            num_pis = vector_piL(piL):step_pi:1;
            M2_ST_pi_m = zeros(numel(num_pis), m_breaks);
            for m = 1:m_breaks
                for ipi = 1:numel(num_pis)
                    S_m = S_Dates_m_M2_Heter(Y, D, Z, X, Mx, beta_H0, datevec, T, num_pis(ipi), m, T_C);
                    S_index = zeros(T, 1);
                    for irow = 1:size(S_m, 1), S_index(S_m(irow, 1):S_m(irow, 2)) = 1; end
                    CT_Z = zeros(size(Z)); CT_Z(S_index == 1, :) = Z(S_index == 1, :);
                    st = legacy_stats(Y, D, Mx * CT_Z, X, T, q, p, a, b, Lag);
                    M2_ST_pi_m(ipi, m) = st.M2;
                end
            end
            [Val_m, ind_m] = max(M2_ST_pi_m, [], 1);
            [Val_m_pi, ind_m_pi] = max(Val_m);
            Sup_M2_piL(piL) = Val_m_pi; Sup_m_piL(piL) = ind_m_pi; Sup_pi_piL(piL) = num_pis(ind_m(ind_m_pi));
            if verbose, fprintf('  pi_L = %.1f: sup M2 = %.4f at m = %d, pi = %.2f\n', vector_piL(piL), Val_m_pi, ind_m_pi, num_pis(ind_m(ind_m_pi))); end
        end
        [~, idx_Sup] = max(Sup_M2_piL);
        pi_hat = Sup_pi_piL(idx_Sup); m_hat = Sup_m_piL(idx_Sup);
        S_hat = S_Dates_m_M2_Heter(Y, D, Z, X, Mx, beta_H0, datevec, T, pi_hat, m_hat, T_C);

        % --- ST statistics on S_hat (legacy lines 259-302) ---
        Index_S = zeros(T, 1);
        for irow = 1:size(S_hat, 1), Index_S(S_hat(irow, 1):S_hat(irow, 2)) = 1; end
        CT_Z = zeros(size(Z)); CT_Z(Index_S == 1, :) = Z(Index_S == 1, :);
        st_ST = legacy_stats(Y, D, Mx * CT_Z, X, T, q, p, a, b, Lag);
        % --- AMS statistics, full sample (legacy lines 318-349) ---
        st_AMS = legacy_stats(Y_orig, D_orig, Mx * Z_orig, X, T, q, p, a, b, Lag);

        stats = [st_ST.AR, st_ST.LM, st_ST.LR, st_ST.AR_hom, st_ST.LM_hom, st_ST.LR_hom, ...
                 st_AMS.AR, st_AMS.LM, st_AMS.LR, st_AMS.AR_hom, st_AMS.LM_hom, st_AMS.LR_hom];

        % --- Critical values and rejections (legacy lines 305-359): rng(11) at the top of the script, then
        %     LR_critical_value(M2, q, alpha) in the order ST, ST_hom, AMS, AMS_hom (1000 draws each). ---
        rng(11);
        LR_draws = {lr_draws(st_ST.M2, q), lr_draws(st_ST.M2_hom, q), lr_draws(st_AMS.M2, q), lr_draws(st_AMS.M2_hom, q)};
        n_a = numel(alpha_grid);
        cvs = NaN(n_a, 12); rej = false(n_a, 12);
        for ia = 1:n_a
            al = alpha_grid(ia);
            cv_AR = chi2inv(1 - al, q); cv_LM = chi2inv(1 - al, 1);
            cv_LR = cellfun(@(d) quantile(d, 1 - al), LR_draws);     % = LR_critical_value with the same draws
            cvs(ia, :) = [cv_AR, cv_LM, cv_LR(1), cv_AR, cv_LM, cv_LR(2), cv_AR, cv_LM, cv_LR(3), cv_AR, cv_LM, cv_LR(4)];
            rej(ia, :) = stats > cvs(ia, :);
        end
        pval = [1 - chi2cdf(st_ST.AR*q, q), 1 - chi2cdf(st_ST.LM, 1), mean(LR_draws{1} >= st_ST.LR), ...
                1 - chi2cdf(st_ST.AR_hom*q, q), 1 - chi2cdf(st_ST.LM_hom, 1), mean(LR_draws{2} >= st_ST.LR_hom), ...
                1 - chi2cdf(st_AMS.AR*q, q), 1 - chi2cdf(st_AMS.LM, 1), mean(LR_draws{3} >= st_AMS.LR), ...
                1 - chi2cdf(st_AMS.AR_hom*q, q), 1 - chi2cdf(st_AMS.LM_hom, 1), mean(LR_draws{4} >= st_AMS.LR_hom)];

        Results.(lbl) = struct('Y', Y_name, 'D', D_name, 'T', T, 'T_C', T_C, 'T_P', T_P, 'datevec', datevec, ...
                               'pi_hat', pi_hat, 'm_hat', m_hat, 'S_hat', S_hat, 'n_S', sum(Index_S), 'Sup_M2_piL', Sup_M2_piL, ...
                               'labels', {labels}, 'stats', stats, 'M2', [st_ST.M2, st_ST.M2_hom, st_AMS.M2, st_AMS.M2_hom], ...
                               'cvs', cvs, 'rej', rej, 'pval', pval);
        if verbose
            fprintf('  S_hat: pi_hat = %.2f, m_hat = %d, |S| = %d of %d (control block [1,%d], treatment block [%d,%d])\n', pi_hat, m_hat, sum(Index_S), T, T_C, T_C+1, T);
            fprintf('  segments:'); fprintf(' [%d,%d]', S_hat'); fprintf('\n');
            fprintf('  %-12s %10s %10s   %s\n', 'statistic', 'value', 'p-value', strjoin(arrayfun(@(a) sprintf('rej@%.2f', a), alpha_grid, 'UniformOutput', false), ' '));
            for k = 1:12
                fprintf('  %-12s %10.4f %10.4f   %s\n', labels{k}, stats(k), pval(k), strjoin(arrayfun(@(r) sprintf('%7d', r), rej(:, k)', 'UniformOutput', false), ' '));
            end
        end
    end

    % --- Export: layout of the published Table 3 ---
    txt = fullfile(folder_name, sprintf('Table3_IdentificationRobust_%s.txt', timestamp_file));
    fid = fopen(txt, 'w');
    fprintf(fid, 'Table 3: Identification-Robust Inference on beta (H0: beta = 0; 1 = rejection, 0 = non-rejection)\n');
    fprintf(fid, 'Legacy script Run_Identi_Robus_Tests_NS18_Heter.m reproduced (see runSection7_Table3.m). T=%d, T_C=%d, T_P=%d.\n\n', T, T_C, T_P);
    pl = fieldnames(Results); pl = pl(startsWith(pl, 'Y_'));
    fprintf(fid, '%-14s', 'alpha'); for ip = 1:numel(pl), fprintf(fid, ' | %-22s', strrep(pl{ip}, 'Y_', '')); end; fprintf(fid, '\n');
    fprintf(fid, '%-14s', ''); for ip = 1:numel(pl), fprintf(fid, ' | %s', strjoin(arrayfun(@(a) sprintf('%6.2f', a), alpha_grid, 'UniformOutput', false), ' ')); end; fprintf(fid, '\n');
    rows = {'LM_T', 8; 'CLR_T', 9; 'LM_T(S_hat)', 2; 'CLR_T(S_hat)', 3; 'AR_T', 7; 'AR_T(S_hat)', 1};
    for r = 1:size(rows, 1)
        fprintf(fid, '%-14s', rows{r, 1});
        for ip = 1:numel(pl)
            fprintf(fid, ' | %s', strjoin(arrayfun(@(x) sprintf('%6d', x), double(Results.(pl{ip}).rej(:, rows{r, 2}))', 'UniformOutput', false), ' '));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\nStatistics, p-values (LR: simulated conditional distribution, 1000 draws, rng(11)) and S_hat:\n');
    for ip = 1:numel(pl)
        R = Results.(pl{ip});
        fprintf(fid, '\n%s: pi_hat = %.2f, m_hat = %d, |S_hat| = %d; segments', pl{ip}, R.pi_hat, R.m_hat, R.n_S); fprintf(fid, ' [%d,%d]', R.S_hat'); fprintf(fid, '\n');
        for k = 1:12, fprintf(fid, '  %-12s %10.4f  p = %.4f  cv(%s) = %s\n', labels{k}, R.stats(k), R.pval(k), strjoin(arrayfun(@(a) sprintf('%.2f', a), alpha_grid, 'UniformOutput', false), '/'), strjoin(arrayfun(@(c) sprintf('%.4f', c), R.cvs(:, k)', 'UniformOutput', false), '/')); end
    end
    fclose(fid);
    save(fullfile(folder_name, sprintf('Table3_IdentificationRobust_%s.mat', timestamp_file)), 'Results');
    if verbose, fprintf('\nSaved: %s\n', txt); end
end

function st = legacy_stats(Y, D, Zbar, X, T, q, p, a, b, Lag)
% Statistics of the legacy script (lines 199-218, 269-302, 318-349): NW HAC covariances with the given lag,
% robust and homoskedastic versions; AR = M1/q.
    y = [Y, D];
    if isempty(X)
        v_hat = y - Zbar * inv(Zbar' * Zbar) * Zbar' * y;
    else
        v_hat = y - Zbar * inv(Zbar' * Zbar) * Zbar' * y - X * inv(X' * X) * X' * y;
    end
    Sigma_v_hat = v_hat' * v_hat / (T - q - p);
    Sigma_N1_hat     = NeweyWest(v_hat * b, Zbar, Lag, 0);
    Sigma_N2star_hat = NeweyWest(v_hat * inv(Sigma_v_hat) * a, Zbar, Lag, 0);
    Sigma_N1N2_hat   = NeweyWest_N1N2((v_hat * b .* Zbar), (v_hat * inv(Sigma_v_hat) * a .* Zbar), Lag, 0);
    Sigma_N2_hat     = Sigma_N2star_hat - Sigma_N1N2_hat * inv(Sigma_N1_hat) * Sigma_N1N2_hat';
    N1 = Sigma_N1_hat^(-1/2) * T^(-1/2) * Zbar' * y * b;
    N2 = Sigma_N2_hat^(-1/2) * (T^(-1/2) * Zbar' * y * inv(Sigma_v_hat) * a - Sigma_N1N2_hat * Sigma_N1_hat^(-1/2) * N1);
    N1_hom = (Zbar' * Zbar)^(-1/2) * Zbar' * y * b * (b' * Sigma_v_hat * b)^(-1/2);
    N2_hom = (Zbar' * Zbar)^(-1/2) * Zbar' * y * inv(Sigma_v_hat) * a * (a' * inv(Sigma_v_hat) * a)^(-1/2);
    [st.AR, st.LM, st.LR, st.M2]                 = ams_from_N(N1, N2, q);
    [st.AR_hom, st.LM_hom, st.LR_hom, st.M2_hom] = ams_from_N(N1_hom, N2_hom, q);
end

function [AR, LM, LR, M2] = ams_from_N(N1, N2, q)
    M1 = N1' * N1; M12 = N1' * N2; M2 = N2' * N2;
    AR = M1 / q; LM = M12^2 / M2;
    LR = 0.5 * (M1 - M2 + sqrt((M1 - M2)^2 + 4 * M12^2));
end

function LR = lr_draws(M2_ST, q)
% The 1000 draws of LR_critical_value.m (legacy_IRI), in its RNG order; the CV at level alpha is their
% (1 - alpha) quantile, as in the legacy function.
    n_reps = 1000; LR = zeros(n_reps, 1);
    for i_reps = 1:n_reps
        M1 = chi2rnd(q);
        M12 = sqrt(M2_ST) * randn;
        LR(i_reps, 1) = 0.5 * (M1 - M2_ST + sqrt((M1 - M2_ST)^2 + 4 * M12^2));
    end
end
