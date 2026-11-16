function Results = runSection7_Table2(h_rule, external_dir, verbose)
% RUNSECTION7_TABLE2  Table 2 (estimation of beta: TSLS and GMM, full sample and pi-hat sample) and Figure 3
% (strongly-identified sub-sample S_hat_{T,OLS}) of the paper, Section 7.2.
%
% Legacy script: pi-LATE\Run_Estimation_IV_GMM_NS18_Heter.m (24/8/2025; README: lines 306 and 502, Figure 3 at
% lines 638-683). Same data, sample and instrument as runSection7_Table1/3 (observations 401:end, control days
% first, D demeaned within group, Z = +-(T/T_.) D).
%   * S_hat_{T,OLS}: from D = 1-day 2-Year yields (DNY2) only, dating of the control sample (Procedure_piL_Full_Sample,
%     m = 5, eps1 = 0.05), selection by the first-stage F (S_Dates_m_Heter, F_ST_Stat with Hac_option = 0) over
%     pi_L in {0.8, 0.9, 1}, pi in pi_L:0.05:1, m = 1..5; pi_hat = pi of the sup. The same S_hat is used for the
%     30-minute policy news. h_rule: the script in the package has h = round(eps1*T_C) = 38 ('T_C'); the published
%     Table 2 and the description of S_hat in the text (excluded: early 2005 and the financial-crisis regime) are
%     reproduced with h = round(eps1*T) = 42 ('T', default), as for Table 1. Both give pi_hat = 0.8, m_hat = 4.
%   * TSLS rows: TSLS_Robust (legacy_mbreak, verbatim; HAC = 0, heteroskedasticity-robust SE with the small-sample
%     factor n/(n-k-1)); standard 95% CI = beta +- 1.96 SE; stars at 10/5/1% (two-sided normal).
%   * GMM rows (Lewis 2022: GMM on the second moments of the two regimes, robust CI by inversion of the subset
%     K-test): computed with the function estimate(eta11, eta12, eta21, eta22) of Lewis's replication code, which is
%     NOT part of this repository. Copy the folder 'estimation' of that package into external_dir (default
%     Section7\external\Lewis2022, to be kept out of git) or pass its path; if estimate.m is not found the GMM rows
%     are left empty and the table is produced with the TSLS rows only. Lewis's code needs the Symbolic Math Toolbox.
%     The legacy call is estimate(Y_C, Y_P, D_C, D_P) (first variable: forward rates; second: policy variable) and
%     the paper reports its first row (coefficient H_12). Both rows of the returned cell are stored.
%   * Figure 3: 1-day changes in 2-Year yields in the control sample and in the policy sample, with the regimes of
%     S_hat (control part) and pi_hat; saved as .png/.eps in the results folder (plotting is guarded: a failure to
%     draw does not stop the numerical part).
%   * Omitted from the script (no effect on the table): Procedure_pi_LATE(...) (printing only), OLS of Y on the
%     fitted D (Estimation_pi_Late / OLS_Robust_full, reported there as a check), the Nakamura-Steinsson covariance
%     ratio, the 3/5/10-year forwards, the loads of nak_data.mat / results_nak.mat / time.mat (unused by estimate).
%
%   Results = runSection7_Table2()                 % h = round(eps1*T), GMM if Section7\external\Lewis2022 exists
%   Results = runSection7_Table2('T_C', 'C:\...\estimation')   % h of the package script, Lewis's code elsewhere
% =========================================================================
    if nargin < 1 || isempty(h_rule), h_rule = 'T'; end
    if nargin < 2 || isempty(external_dir), external_dir = fullfile('external', 'Lewis2022'); end
    if nargin < 3 || isempty(verbose), verbose = true; end
    restoredefaultpath;
    addpath('..\m_break'); addpath('legacy_mbreak');
    have_gmm = false;
    if exist(external_dir, 'dir') == 7
        addpath(external_dir);
        have_gmm = exist('estimate', 'file') == 2 && exist('InvertTest', 'file') == 2 && exist('GMM', 'file') == 2;
    end
    warning('off', 'all');

    folder_name = 'Paper_Exhibits_Results_SEC7_TABLE2_FIG3';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

    dta = readtable(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv'));
    FOMC = dta.FOMCused(401:end); T = numel(FOMC); T_C = sum(FOMC == 0); T_P = sum(FOMC == 1);
    dates_str = read_first_column(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv'));   % date_daily as text
    dates_str = dates_str(401:end);
    dates_ord = [dates_str(FOMC == 0); dates_str(FOMC == 1)];           % re-ordered sample [control; treatment]
    DNY2 = dta.DNY2(401:end); DNF2 = dta.DNF2(401:end); DRF2 = dta.DRF2(401:end); it = dta.path_intra_wide(401:end);
    assert(~any(isnan([DNY2; DNF2; DRF2; it])), 'NaN in the sample 401:end.');
    reord = @(v) [v(FOMC == 0); v(FOMC == 1)];
    dem   = @(v) [v(FOMC == 0) - mean(v(FOMC == 0)); v(FOMC == 1) - mean(v(FOMC == 1))];
    instr = @(v) [-T/T_C * (v(FOMC == 0) - mean(v(FOMC == 0))); T/T_P * (v(FOMC == 1) - mean(v(FOMC == 1)))];

    % --- S_hat_{T,OLS} from the 1-day yields (legacy lines 494-624) ---
    y = dem(DNY2); z = instr(DNY2); x = []; q = 1; p = 0; m_breaks = 5; eps1 = 0.05;
    if strcmpi(h_rule, 'T_C'), h = round(eps1 * T_C); else, h = round(eps1 * T); end
    [~, datevec, ~] = Procedure_piL_Full_Sample(T_C, y(1:T_C), z(1:T_C), q, m_breaks, h, eps1, 1, 0, 1, 1, 1, 0, x, p, 0.0001, 20, 0, 0, 1, 1);
    vector_piL = 0.8:0.1:1; step_pi = 0.05; max_breaks = 5; Hac_option = 0;
    nL = numel(vector_piL); Sup_F_piL = NaN(nL, 1); Sup_m_piL = NaN(nL, 1); Sup_pi_piL = NaN(nL, 1);
    for iL = 1:nL
        num_pis = vector_piL(iL):step_pi:1; F_ST = zeros(numel(num_pis), max_breaks);
        for m = 1:max_breaks
            for ipi = 1:numel(num_pis)
                Tpi = floor(num_pis(ipi) * T);
                if vector_piL(iL) == 1, S_m = [1, T]; else, S_m = S_Dates_m_Heter(y, z, x, datevec, T, num_pis(ipi), m, T_C); end
                S_index = zeros(T, 1); for irow = 1:size(S_m, 1), S_index(S_m(irow, 1):S_m(irow, 2)) = 1; end
                F_ST(ipi, m) = F_ST_Stat(y(S_index == 1), z(S_index == 1), Tpi, q, p, Hac_option);
            end
        end
        [Val_m, ind_m] = max(F_ST, [], 1); [Sup_F_piL(iL), Sup_m_piL(iL)] = max(Val_m); Sup_pi_piL(iL) = num_pis(ind_m(Sup_m_piL(iL)));
    end
    [~, idx] = max(Sup_F_piL);
    pi_hat = Sup_pi_piL(idx); m_hat = Sup_m_piL(idx);
    if vector_piL(idx) == 1, S_hat = [1, T]; else, S_hat = S_Dates_m_Heter(y, z, x, datevec, T, pi_hat, m_hat, T_C); end
    Index_S = zeros(T, 1); for irow = 1:size(S_hat, 1), Index_S(S_hat(irow, 1):S_hat(irow, 2)) = 1; end
    ctrl_rows = S_hat(S_hat(:, 2) <= T_C, :);
    if verbose
        fprintf('S_hat_{T,OLS} (D = DNY2, h = %d): pi_hat = %.2f, m_hat = %d, |S| = %d of %d; sup F per pi_L %s = %s\n', h, pi_hat, m_hat, sum(Index_S), T, mat2str(vector_piL), mat2str(round(Sup_F_piL' * 100) / 100));
        fprintf('  control regimes in S_hat (dates):'); for r = 1:size(ctrl_rows, 1), fprintf(' [%s, %s]', dates_ord{ctrl_rows(r, 1)}, dates_ord{ctrl_rows(r, 2)}); end; fprintf('\n');
        excl = find(Index_S(1:T_C) == 0); if ~isempty(excl), fprintf('  excluded control days: %d (from %s to %s)\n', numel(excl), dates_ord{excl(1)}, dates_ord{excl(end)}); end
    end

    % --- TSLS and GMM, full sample and S_hat (legacy lines 696-760, 760-811, 812-880) ---
    Dvars = {'it', '30-minute policy news', it; 'DNY2', '1-day 2-Year yields', DNY2};
    Yvars = {'DNF2', 'Nominal', DNF2; 'DRF2', 'Real', DRF2};
    Results = struct('T', T, 'T_C', T_C, 'T_P', T_P, 'h', h, 'h_rule', h_rule, 'pi_hat', pi_hat, 'm_hat', m_hat, 'S_hat', S_hat, ...
                     'Index_S', Index_S, 'Sup_F_piL', Sup_F_piL, 'vector_piL', vector_piL, 'have_gmm', have_gmm, 'external_dir', external_dir);
    for iD = 1:2
        D_raw = Dvars{iD, 3}; D_dem = dem(D_raw); Z_dem = instr(D_raw);
        for iY = 1:2
            Y_raw = Yvars{iY, 3}; Y_ord = reord(Y_raw);
            lbl = sprintf('%s_%s', Yvars{iY, 1}, Dvars{iD, 1});
            R = struct('Y', Yvars{iY, 1}, 'D', Dvars{iD, 1});
            for smp = {'full', 'Shat'}
                if strcmp(smp{1}, 'full'), sel = true(T, 1); else, sel = Index_S == 1; end
                [b, se, pv] = TSLS_Robust(Y_ord(sel), D_dem(sel), Z_dem(sel), Hac_option);
                R.(['tsls_' smp{1}]) = struct('beta', b, 'se', se, 'pval', pv, 'ci95', [b - 1.96*se, b + 1.96*se], 'stars', stars(b / se));
                if have_gmm
                    cC = sel(1:T_C); cP = sel(T_C+1:end);
                    Y_C = Y_raw(FOMC == 0); Y_P = Y_raw(FOMC == 1); D_C = D_raw(FOMC == 0); D_P = D_raw(FOMC == 1);
                    try
                        Tbl = estimate(Y_C(cC), Y_P(cP), D_C(cC), D_P(cP));
                        R.(['gmm_' smp{1}]) = struct('table', {Tbl}, 'beta', Tbl{1, 1}, 'ci95_std', Tbl{1, 2}, 'ci95_robust', Tbl{1, 8}, 'tstat', Tbl{1, 7});
                    catch ME
                        R.(['gmm_' smp{1}]) = struct('error', ME.message);
                    end
                end
            end
            Results.(lbl) = R;
        end
    end

    % --- Export (layout of the published Table 2) ---
    txt = fullfile(folder_name, sprintf('Table2_Estimation_h%s_%s.txt', h_rule, timestamp_file));
    fid = fopen(txt, 'w');
    fprintf(fid, 'Table 2: Estimation of beta (legacy Run_Estimation_IV_GMM_NS18_Heter.m reproduced; see runSection7_Table2.m)\n');
    fprintf(fid, 'T=%d, T_C=%d, T_P=%d; S_hat_{T,OLS} from D = DNY2 with h = round(eps1*%s) = %d: pi_hat = %.2f, m_hat = %d, |S| = %d\n', T, T_C, T_P, h_rule, h, pi_hat, m_hat, sum(Index_S));
    fprintf(fid, 'S_hat regimes (re-ordered sample, control = 1..%d):', T_C); fprintf(fid, ' [%d,%d]', S_hat'); fprintf(fid, '\n');
    fprintf(fid, 'Control regimes (dates):'); for r = 1:size(ctrl_rows, 1), fprintf(fid, ' [%s, %s]', dates_ord{ctrl_rows(r, 1)}, dates_ord{ctrl_rows(r, 2)}); end; fprintf(fid, '\n');
    if ~have_gmm, fprintf(fid, 'GMM rows: not computed (Lewis (2022) code not found in %s).\n', external_dir); end
    cols = {'DNF2_it', 'DRF2_it', 'DNF2_DNY2', 'DRF2_DNY2'};
    fprintf(fid, '\n%-22s %-26s %-26s %-26s %-26s\n', '', '30-min Nominal', '30-min Real', '1-day Nominal', '1-day Real');
    for smp = {'full', 'Shat'}
        if strcmp(smp{1}, 'full'), fprintf(fid, '--- Full Sample ---\n'); else, fprintf(fid, '--- pi-sample S_hat_{T,OLS}, pi_hat = %.2f ---\n', pi_hat); end
        fprintf(fid, '%-22s', 'TSLS beta'); for c = 1:4, t = Results.(cols{c}).(['tsls_' smp{1}]); fprintf(fid, ' %-26s', sprintf('%.2f%s', t.beta, t.stars)); end; fprintf(fid, '\n');
        fprintf(fid, '%-22s', '   standard CI'); for c = 1:4, t = Results.(cols{c}).(['tsls_' smp{1}]); fprintf(fid, ' %-26s', sprintf('[%.2f, %.2f]', t.ci95(1), t.ci95(2))); end; fprintf(fid, '\n');
        if have_gmm
            fprintf(fid, '%-22s', 'GMM beta'); for c = 1:4, g = Results.(cols{c}).(['gmm_' smp{1}]); if isfield(g, 'beta'), fprintf(fid, ' %-26s', sprintf('%.2f%s', g.beta, stars(g.tstat))); else, fprintf(fid, ' %-26s', 'error'); end; end; fprintf(fid, '\n');
            fprintf(fid, '%-22s', '   standard CI'); for c = 1:4, g = Results.(cols{c}).(['gmm_' smp{1}]); if isfield(g, 'ci95_std'), fprintf(fid, ' %-26s', g.ci95_std); else, fprintf(fid, ' %-26s', '-'); end; end; fprintf(fid, '\n');
            fprintf(fid, '%-22s', '   robust CI (K-test)'); for c = 1:4, g = Results.(cols{c}).(['gmm_' smp{1}]); if isfield(g, 'ci95_robust'), fprintf(fid, ' %-26s', g.ci95_robust); else, fprintf(fid, ' %-26s', '-'); end; end; fprintf(fid, '\n');
        end
    end
    fprintf(fid, '\nPublished TSLS: full 1.10** [0.17,2.02] | 0.96*** [0.41,1.51] | 1.14*** [0.83,1.45] | 0.97*** [0.40,1.565]; pi-sample 1.11** [0.19,2.02] | 0.97*** [0.42,1.51] | 1.13*** [0.92,1.30] | 0.92*** [0.56,1.28]\n');
    fprintf(fid, 'Published GMM:  full 1.07** [0.17,1.98] rob [0.27,3.25] | 0.94*** [0.36,1.51] rob [0.44,2.38] | -0.27 [-4.90,4.36] rob [-77.27,0.94] | 1.31 [-3.74,6.35] rob [-253.70,1.92]; pi-sample 1.07** [0.17,1.96] | 0.94*** [0.38,1.50] | 0.65* [-0.02,1.31] | 0.86** [0.29,1.43]\n');
    fclose(fid);
    save(fullfile(folder_name, sprintf('Table2_Estimation_h%s_%s.mat', h_rule, timestamp_file)), 'Results');
    if verbose, fprintf('Saved: %s\n', txt); end

    % --- Figure 3 (legacy lines 638-683): yields in the control and policy samples, S_hat regimes shaded ---
    try
        D_ord = reord(DNY2);
        fh = figure('Visible', 'off');
        subplot(2, 1, 1); hold on;
        for r = 1:size(ctrl_rows, 1)
            patch([ctrl_rows(r,1) ctrl_rows(r,2) ctrl_rows(r,2) ctrl_rows(r,1)], [min(D_ord(1:T_C)) min(D_ord(1:T_C)) max(D_ord(1:T_C)) max(D_ord(1:T_C))], [1 0.85 0.85], 'EdgeColor', [0.8 0 0]);
        end
        plot(1:T_C, D_ord(1:T_C), 'b-'); xlim([1 T_C]);
        title(sprintf('Control sample: 1-day changes in 2-Year yields; shaded = S_hat_{T,OLS}, pi_hat = %.2f', pi_hat));
        tk = round(linspace(1, T_C, 6)); set(gca, 'XTick', tk, 'XTickLabel', dates_ord(tk));
        subplot(2, 1, 2); plot(1:T_P, D_ord(T_C+1:end), 'b-'); xlim([1 T_P]);
        title('Policy sample (FOMC days): 1-day changes in 2-Year yields (all in S_hat)');
        tk = round(linspace(1, T_P, 6)); set(gca, 'XTick', tk, 'XTickLabel', dates_ord(T_C + tk));
        print(fh, fullfile(folder_name, sprintf('Figure3_Shat_h%s_%s.png', h_rule, timestamp_file)), '-dpng');
        close(fh);
    catch ME
        if verbose, fprintf('Figure 3 not drawn: %s\n', ME.message); end
    end
end

function s = stars(tval)
    a = abs(tval);
    if a > 2.5758, s = '***'; elseif a > 1.96, s = '**'; elseif a > 1.6449, s = '*'; else, s = ''; end
end

function c = read_first_column(fname)
% First column of a csv file as a cell array of strings (header skipped); independent of readtable's type choices.
    fid = fopen(fname, 'r'); fgetl(fid); c = {};
    while true
        l = fgetl(fid); if ~ischar(l), break; end
        k = find(l == ',', 1); if isempty(k), c{end+1, 1} = strtrim(l); else, c{end+1, 1} = strtrim(l(1:k-1)); end
    end
    fclose(fid);
end
