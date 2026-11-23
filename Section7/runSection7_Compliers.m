function Results = runSection7_Compliers(lag_rule, S_hat_fig, verbose)
% RUNSECTION7_COMPLIERS  Section 7.4: identification of compliers (Figure 4, share of compliers) and test of the
% exclusion restriction (t_exclusion for the real and nominal 2-Year forward rates).
%
% Legacy scripts (pi-LATE\Compliers\): Run_Compliers_Heter_2sided_Figure_D2.m (Figure 4 and share of compliers;
% README) and Run_Test_Exclusion_Restriction_2Y{R,N}F_HAC_nonCompliers_YD_D2.m (exclusion restriction, README
% "lines 425-430 / 451-456"). The three scripts share the classification code; they differ in the Newey-West lag
% of the compliers t-test: floor(n0^(1/3)) = 4 in the Figure script, 0 in the exclusion scripts. Both are available
% here (lag_rule 'n0' or 0) and both exclusion statistics are reported; the published t_exclusion = 0.01 (real)
% and 0.60 (nominal) are those of the exclusion scripts.
%
% DESIGN (Theorem 2.1 with heteroskedasticity-based identification):
%   * Data: ..\m_break\dataForMatlabAltEstTick.csv, observations 401:end, D_t = 1-day change in 2-Year yields (DNY2);
%     the scripts drop the last 13 observations (T = 823, T_C = 749, T_P = 74). Original time order (not re-ordered).
%   * Rolling two-sided windows of squared D: control days, n0 = 2*50 + 1 = 101 control observations around t
%     (sigma2_C,t,n0; computed only when 50 control days exist on each side); FOMC days, n1 = 2*7 + 1 = 15 FOMC
%     observations around t (sigma2_P,t,n1); for a control day t the policy window is taken at the closest FOMC day
%     s*(t). Differences: diff_P(t) = sigma2_P,t,n1 - sigma2_C,t-1,n0 (t in P), diff_C(t) = sigma2_P,s*(t),n1 - sigma2_C,t,n0 (t in C).
%   * Compliers (10% one-sided): t_t = sqrt(n0) diff / sqrt(J_HAC (1 + n0/n1)) > 1.28, with J_HAC = Newey-West of
%     the demeaned squared D over the control window (lag per lag_rule), window at t-1 for t in P and at t for t in C.
%     The script also classifies by diff > mean(diff) (reported as 'point' classification).
%   * Exclusion restriction: Y_t = D_t * Y_t (forward rate change times D) on the non-compliers; t_exclusion =
%     |sqrt(|NC_C|) (mean_C - mean_P)| / sqrt(J_HAC,NC_C (1 + |NC_C|/|NC_P|)), J_HAC with lag floor(|NC_C|^(1/3));
%     rejection at 2.54 (one-sided 1%... as in the script). Also the always-/never-takers split of the scripts
%     (non-compliers with D^2 above/below their mean) with the same statistic.
%   * Figure 4: D in the control and policy samples, markers for compliers (green) / non-compliers (red), rectangles
%     for S_hat_{T,OLS}. S_hat_fig: regimes in the re-ordered sample; default = the one hard-coded in the exclusion
%     scripts ([1 98; 145 302; 425 472; 473 762; 763 836], pi_hat = 0.8, T_C = 762), i.e. the published S_hat, which
%     is also what runSection7_Table2 returns with its default h_rule 'T_C' (h = round(eps1*T_C) = 38).
%
%   Results = runSection7_Compliers()           % lag floor(n0^(1/3)) = 4 for the compliers t-test (Figure script)
%   Results = runSection7_Compliers(0)          % lag 0 (exclusion scripts)
% =========================================================================
    if nargin < 1 || isempty(lag_rule), lag_rule = 'n0'; end
    if nargin < 2 || isempty(S_hat_fig), S_hat_fig = [1 98; 145 302; 425 472; 473 762; 763 836]; end
    if nargin < 3 || isempty(verbose), verbose = true; end
    restoredefaultpath; addpath('..\m_break'); warning('off', 'all');
    folder_name = 'Paper_Exhibits_Results_SEC7_COMPLIERS_FIG4';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

    dta = readtable(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv'));
    dates = read_first_column(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv')); dates = dates(401:end);
    FOMC = dta.FOMCused(401:end); DNY2 = dta.DNY2(401:end); DNF2 = dta.DNF2(401:end); DRF2 = dta.DRF2(401:end);
    % legacy: last 13 observations dropped
    FOMC = FOMC(1:end-13); D = DNY2(1:end-13); dates = dates(1:end-13);
    T = numel(FOMC); T_C = sum(FOMC == 0); T_P = sum(FOMC == 1);
    n0_os = 50; n0 = 2*n0_os + 1; n1_os = 7; n1 = 2*n1_os + 1;
    if ischar(lag_rule), Lag = floor(n0^(1/3)); else, Lag = lag_rule; end

    % --- rolling window variances (legacy loops, verbatim logic) ---
    D_C_n0 = NaN(T, 1); D_P_n1 = NaN(T, 1); D_P_star_n1 = NaN(T, 1);
    W_C = NaN(n0, T);
    for t = 1:T
        if FOMC(t) == 0
            nPb = sum(FOMC(1:t-1) == 1); nPa = sum(FOMC(t+1:end) == 1);
            if t - n0_os - nPb - 1 >= 0 && t + n0_os + nPa <= T
                idx = find(FOMC(1:t) == 0); Dn = D(idx); w = Dn(end-n0_os:end);
                idx = t + find(FOMC(t+1:end) == 0); Dn = D(idx); w = [w; Dn(1:n0_os)];
                W_C(:, t) = w; D_C_n0(t) = mean(w.^2);
            end
        end
    end
    for t = 1:T
        if FOMC(t) == 1
            nCb = sum(FOMC(1:t-1) == 0); nCa = sum(FOMC(t+1:end) == 0);
            if t - n1_os - nCb - 1 >= 0 && t + n1_os + nCa <= T
                idx = find(FOMC(1:t) == 1); Dn = D(idx); w = Dn(end-n1_os:end);
                idx = t + find(FOMC(t+1:end) == 1); Dn = D(idx); w = [w; Dn(1:n1_os)];
                D_P_n1(t) = mean(w.^2);
            end
        end
    end
    for t = 1:T
        if FOMC(t) == 0
            idxP = find(FOMC == 1); [~, k] = min(abs(idxP - t)); s_star = idxP(k);
            nCb = sum(FOMC(1:s_star-1) == 0); nCa = sum(FOMC(s_star+1:end) == 0);
            if s_star - n1_os - nCb - 1 >= 0 && s_star + n1_os + nCa <= T
                idx = find(FOMC(1:s_star) == 1); Dn = D(idx); w = Dn(end-n1_os:end);
                idx = s_star + find(FOMC(s_star+1:end) == 1); Dn = D(idx); w = [w; Dn(1:n1_os)];
                D_P_star_n1(t) = mean(w.^2);
            end
        end
    end
    diff_P = NaN(T, 1); diff_P(2:T) = D_P_n1(2:end) - D_C_n0(1:end-1);
    diff_C = D_P_star_n1 - D_C_n0;
    diff_full = diff_C; diff_full(FOMC == 1) = diff_P(FOMC == 1);
    c = mean(diff_full, 'omitnan');                      % legacy: mean(diff_full) (NaN if any NaN; used only for the point classification)
    comp_point_P = double(diff_P > c); comp_point_P(isnan(diff_P)) = NaN;
    comp_point_C = double(diff_C > c); comp_point_C(isnan(diff_C)) = NaN;

    % --- t-tests (10% one-sided, cv 1.28) ---
    t_P = NaN(T, 1); t_C = NaN(T, 1); rej_P = NaN(T, 1); rej_C = NaN(T, 1);
    for t = 1:T
        if ~isnan(diff_P(t))
            r = W_C(:, t-1).^2 - mean(W_C(:, t-1).^2); J = NeweyWest(r, ones(n0, 1), Lag, 0);
            t_P(t) = sqrt(n0) * diff_P(t) / sqrt(J * (1 + n0/n1)); rej_P(t) = double(t_P(t) > 1.28);
        end
        if ~isnan(diff_C(t))
            r = W_C(:, t).^2 - mean(W_C(:, t).^2); J = NeweyWest(r, ones(n0, 1), Lag, 0);
            t_C(t) = sqrt(n0) * diff_C(t) / sqrt(J * (1 + n0/n1)); rej_C(t) = double(t_C(t) > 1.28);
        end
    end
    nC = nnz(~isnan(rej_C)); nP = nnz(~isnan(rej_P));
    frac_C = nnz(rej_C == 1) / nC; frac_P = nnz(rej_P == 1) / nP; frac_all = (nnz(rej_C == 1) + nnz(rej_P == 1)) / (nC + nP);
    idx_nonC_C = find(rej_C == 0); idx_nonC_P = find(rej_P == 0);
    comp_C_in_C = rej_C(FOMC == 0); comp_P_in_P = rej_P(FOMC == 1);

    % --- exclusion restriction, for Y = real and nominal 2-Year forwards (legacy lines 393-456) ---
    Ys = {'DRF2', DRF2(1:end-13), 'real'; 'DNF2', DNF2(1:end-13), 'nominal'};
    D2nonC = [D(idx_nonC_C); D(idx_nonC_P)].^2; m_D2 = mean(D2nonC);
    isAT_C = D(idx_nonC_C).^2 - m_D2 >= 0; isAT_P = D(idx_nonC_P).^2 - m_D2 >= 0;
    Excl = struct();
    for iy = 1:2
        Y = Ys{iy, 2};
        Excl.(Ys{iy, 1}) = struct('label', Ys{iy, 3}, ...
            'nonC', excl_stat(Y(idx_nonC_C) .* D(idx_nonC_C), Y(idx_nonC_P) .* D(idx_nonC_P)), ...
            'AT',   excl_stat(Y(idx_nonC_C(isAT_C)) .* D(idx_nonC_C(isAT_C)), Y(idx_nonC_P(isAT_P)) .* D(idx_nonC_P(isAT_P))), ...
            'NT',   excl_stat(Y(idx_nonC_C(~isAT_C)) .* D(idx_nonC_C(~isAT_C)), Y(idx_nonC_P(~isAT_P)) .* D(idx_nonC_P(~isAT_P))));
    end

    Results = struct('T', T, 'T_C', T_C, 'T_P', T_P, 'n0', n0, 'n1', n1, 'Lag', Lag, 'lag_rule', lag_rule, 'dates', {dates}, 'FOMC', FOMC, 'D', D, ...
                     'sigma2_C', D_C_n0, 'sigma2_P', D_P_n1, 'sigma2_P_star', D_P_star_n1, 'diff_P', diff_P, 'diff_C', diff_C, 'c', c, ...
                     'comp_point_P', comp_point_P, 'comp_point_C', comp_point_C, 't_P', t_P, 't_C', t_C, 'rej_P', rej_P, 'rej_C', rej_C, ...
                     'frac_compliers_C', frac_C, 'frac_compliers_P', frac_P, 'frac_compliers_all', frac_all, 'n_classified', nC + nP, ...
                     'n_noncompliers', numel(idx_nonC_C) + numel(idx_nonC_P), 'frac_AT', (nnz(isAT_C) + nnz(isAT_P)) / (numel(idx_nonC_C) + numel(idx_nonC_P)), ...
                     'Exclusion', Excl, 'S_hat_fig', S_hat_fig);
    nonC_days = dates(sort([idx_nonC_C; idx_nonC_P]));
    if verbose
        fprintf('Compliers (lag %d): classified %d of %d (C %d, P %d); share of compliers %.3f (C %.3f, P %.3f); non-compliers %d\n', Lag, nC + nP, T, nC, nP, frac_all, frac_C, frac_P, Results.n_noncompliers);
        fprintf('  non-compliers from %s to %s; by year:', nonC_days{1}, nonC_days{end});
        yrs = cellfun(@(s) s(end-3:end), nonC_days, 'UniformOutput', false); [u, ~, j] = unique(yrs); cnt = accumarray(j, 1);
        for k = 1:numel(u), fprintf(' %s:%d', u{k}, cnt(k)); end; fprintf('\n');
        for iy = 1:2
            E = Excl.(Ys{iy, 1});
            fprintf('  exclusion, %s forwards: t = %.3f (|NC_C| = %d, |NC_P| = %d, lag %d); always-takers t = %.3f; never-takers t = %.3f\n', Ys{iy, 3}, E.nonC.t, E.nonC.nC, E.nonC.nP, E.nonC.lag, E.AT.t, E.NT.t);
        end
    end

    % --- export ---
    txt = fullfile(folder_name, sprintf('Compliers_Exclusion_lag%d_%s.txt', Lag, timestamp_file));
    fid = fopen(txt, 'w');
    fprintf(fid, 'Section 7.4: compliers and exclusion restriction (legacy Compliers\\Run_Compliers_Heter_2sided_Figure_D2.m and Run_Test_Exclusion_Restriction_2Y*F_HAC_nonCompliers_YD_D2.m)\n');
    fprintf(fid, 'T=%d (last 13 obs. dropped), T_C=%d, T_P=%d, n0=%d, n1=%d, Newey-West lag of the compliers t-test = %d, cv 1.28 (10%% one-sided)\n', T, T_C, T_P, n0, n1, Lag);
    fprintf(fid, 'Share of compliers: %.4f overall (%d classified), control %.4f, policy %.4f. Published: about 62%%.\n', frac_all, nC + nP, frac_C, frac_P);
    fprintf(fid, 'Non-compliers: %d (%s to %s)\n', Results.n_noncompliers, nonC_days{1}, nonC_days{end});
    for iy = 1:2
        E = Excl.(Ys{iy, 1});
        fprintf(fid, 'Exclusion restriction, %s 2-Year forwards: t_exclusion = %.4f (|NC_C| = %d, |NC_P| = %d, lag %d, reject at 2.54: %d); always-takers t = %.4f; never-takers t = %.4f\n', Ys{iy, 3}, E.nonC.t, E.nonC.nC, E.nonC.nP, E.nonC.lag, E.nonC.t > 2.54, E.AT.t, E.NT.t);
    end
    fprintf(fid, 'Published: t_exclusion = 0.01 (real), 0.60 (nominal), exclusion scripts (lag 0 in the compliers test).\n');
    fprintf(fid, '\nDate, FOMC, D, diff, t-stat, complier(1)/non-complier(0)/NaN:\n');
    for t = 1:T
        if FOMC(t) == 1, dd = diff_P(t); tt = t_P(t); rr = rej_P(t); else, dd = diff_C(t); tt = t_C(t); rr = rej_C(t); end
        fprintf(fid, '%s %d %9.4f %9.4f %8.3f %g\n', dates{t}, FOMC(t), D(t), dd, tt, rr);
    end
    fclose(fid);
    save(fullfile(folder_name, sprintf('Compliers_Exclusion_lag%d_%s.mat', Lag, timestamp_file)), 'Results');
    if verbose, fprintf('Saved: %s\n', txt); end

    % --- Figure 4 ---
    try
        D_C = D(FOMC == 0); D_P = D(FOMC == 1); idxC = find(FOMC == 0); idxP = find(FOMC == 1);
        fh = figure('Visible', 'off');
        subplot(2, 1, 1); hold on;
        ctrl_rows = S_hat_fig(S_hat_fig(:, 2) <= 762, :);
        for r = 1:size(ctrl_rows, 1)
            a = min(ctrl_rows(r, 1), T_C); b = min(ctrl_rows(r, 2), T_C);
            patch([a b b a], [min(D_C) min(D_C) max(D_C) max(D_C)], [1 0.9 0.8], 'EdgeColor', [0.85 0.33 0.1], 'LineStyle', '--');
        end
        plot(1:T_C, D_C, 'b-');
        cc = comp_C_in_C; plot(find(cc == 1), min(D_C) - 0.05, 'o', 'MarkerFaceColor', [0.47 0.77 0.19], 'MarkerEdgeColor', 'none', 'MarkerSize', 3);
        plot(find(cc == 0), min(D_C) - 0.05, 'o', 'MarkerFaceColor', [0.74 0.08 0.08], 'MarkerEdgeColor', 'none', 'MarkerSize', 3);
        xlim([1 T_C]); title('Control sample: D_t (2-Year yields), S_hat regimes, compliers (green) / non-compliers (red)');
        tk = round(linspace(1, T_C, 6)); set(gca, 'XTick', tk, 'XTickLabel', dates(idxC(tk)));
        subplot(2, 1, 2); hold on; plot(1:T_P, D_P, 'b-');
        cp = comp_P_in_P; plot(find(cp == 1), min(D_C) - 0.05, 'o', 'MarkerFaceColor', [0.47 0.77 0.19], 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
        plot(find(cp == 0), min(D_C) - 0.05, 'o', 'MarkerFaceColor', [0.74 0.08 0.08], 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
        xlim([1 T_P]); title('Policy sample (FOMC days)');
        tk = round(linspace(1, T_P, 6)); set(gca, 'XTick', tk, 'XTickLabel', dates(idxP(tk)));
        print(fh, fullfile(folder_name, sprintf('Figure4_Compliers_lag%d_%s.png', Lag, timestamp_file)), '-dpng'); close(fh);
    catch ME
        if verbose, fprintf('Figure 4 not drawn: %s\n', ME.message); end
    end
end

function E = excl_stat(YC, YP)
% t_exclusion of the legacy scripts: |sqrt(nC) (mean_C - mean_P)| / sqrt(J_HAC(Y_C) (1 + nC/nP)), lag floor(nC^(1/3)).
    nC = numel(YC); nP = numel(YP); lag = floor(nC^(1/3));
    J = NeweyWest(YC - mean(YC), ones(nC, 1), lag, 0);
    E = struct('t', abs(sqrt(nC) * (mean(YC) - mean(YP))) / sqrt(J * (1 + nC/nP)), 'nC', nC, 'nP', nP, 'lag', lag, 'mean_C', mean(YC), 'mean_P', mean(YP));
end

function c = read_first_column(fname)
    fid = fopen(fname, 'r'); fgetl(fid); c = {};
    while true
        l = fgetl(fid); if ~ischar(l), break; end
        k = find(l == ',', 1); if isempty(k), c{end+1, 1} = strtrim(l); else, c{end+1, 1} = strtrim(l(1:k-1)); end
    end
    fclose(fid);
end
