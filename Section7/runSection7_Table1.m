function Results = runSection7_Table1(h_rule, verbose)
% RUNSECTION7_TABLE1  Table 1 of the paper (Section 7.1): F*_T tests for identification failure over all
% sub-populations, for pi_L in {0.6, 0.7, 0.8, 0.9, 1}.
%
% Legacy script: pi-LATE\Run_Fstar_Test_NS18_Heter.m (5/3/2025; README: line 306 for the 30-minute results,
% line 502 for the 1-day results). Same data, sample and instrument as runSection7_Table3 (observations 401:end,
% control days first, D demeaned within group, Z = +-(T/T_.) D). As in the legacy script: pi_L in 0.6:0.1:1,
% pi in pi_L:0.05:1, m = 1..5, selection of S by the first-stage F statistic (S_Dates_m_Heter, m_break) instead
% of M2, F*_T = sup over (pi, m) of F(S) with the heteroskedasticity-only variance (Hac_option = 0:
% OLS_Robust_VarXe with HAC = 0) and the denominator q*(Tpi - p - q), Tpi = floor(pi*T); for pi_L = 1 the
% statistic is the full-sample F. The treatment block is always in S (S_Dates_m_Heter). The legacy call
% Procedure_pi_LATE(...) (printing only, no output used) is omitted.
%
% h_rule: minimal segment length of the Bai-Perron dating on the control sample. The script in the package has
%   h = round(eps1*T_C) = 38 ('T_C'); the published Table 1 is reproduced (1-day column, pi_L = 0.8, 0.9, 1 exactly)
%   with h = round(eps1*T) = 42 ('T', default), the rule of the Table 3 script. Both are available; the full-sample
%   value (pi_L = 1) does not depend on h.
% Critical values: the published Table 4 (q = 1; Section4\cvals\RunMc_CVs_Fstat.m), reported for reference.
% Published 30-minute row: "10^4 x 7.42" etc.; the statistic is 7 416 (same code and scaling as the 1-day
% value 8.09, reproduced exactly), i.e. the published scale factor reads 10^3.
%
%   Results = runSection7_Table1()            % h = round(eps1*T), both D's: 30-minute policy news and 1-day DNY2
%   Results = runSection7_Table1('T_C')       % h = round(eps1*T_C), as in the package script
% =========================================================================
    if nargin < 1 || isempty(h_rule), h_rule = 'T'; end
    if nargin < 2 || isempty(verbose), verbose = true; end
    restoredefaultpath;
    addpath('..\m_break');          % Procedure_piL_Full_Sample, S_Dates_m_Heter, OLS_Robust_VarXe, NeweyWest
    addpath('legacy_mbreak');       % F_ST_Stat (legacy pi-LATE\m_break, verbatim)
    warning('off', 'all');

    folder_name = 'Paper_Exhibits_Results_SEC7_TABLE1';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

    dta = readtable(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv'));
    FOMC = dta.FOMCused(401:end);
    T = numel(FOMC); T_C = sum(FOMC == 0); T_P = sum(FOMC == 1);

    D_list = {'path_intra_wide', '30-minute policy news'; 'DNY2', '1-day nominal Treasury yields'};
    vector_piL = 0.6:0.1:1; step_pi = 0.05; max_breaks = 5; Hac_option = 0;
    cv01 = [17.89, 15.58, 13.11, 10.66, 6.46];     % published Table 4, q = 1, 1%
    cv05 = [12.74, 10.70, 8.88, 7.02, 3.85];       % published Table 4, q = 1, 5%
    % legacy options of the dating procedure
    eps1 = 0.05; q = 1; p = 0; m_breaks = 5; x = [];
    if strcmpi(h_rule, 'T_C'), h = round(eps1 * T_C); else, h = round(eps1 * T); end
    fixb = 0; betaini = 0; maxi = 20; printd = 0; eps = 0.0001;
    robust = 1; prewhit = 0; hetdat = 1; hetvar = 1; hetomega = 1; hetq = 1; doglobal = 1;

    Results = struct('T', T, 'T_C', T_C, 'T_P', T_P, 'h', h, 'h_rule', h_rule, 'vector_piL', vector_piL, 'cv01', cv01, 'cv05', cv05);
    for iD = 1:size(D_list, 1)
        D_name = D_list{iD, 1};
        D_raw = dta.(D_name)(401:end);
        assert(~any(isnan(D_raw)), 'NaN in %s over the sample 401:end.', D_name);
        Dc = D_raw(FOMC == 0); Dp = D_raw(FOMC == 1);
        y = [Dc - mean(Dc, 1); Dp - mean(Dp, 1)];
        z = [-T/T_C * (Dc - mean(Dc, 1)); T/T_P * (Dp - mean(Dp, 1))];
        if verbose, fprintf('\n=== Table 1: D = %s (%s) ===\n', D_name, D_list{iD, 2}); end
        [~, datevec, ~] = Procedure_piL_Full_Sample(T_C, y(1:T_C), z(1:T_C), q, m_breaks, h, eps1, robust, prewhit, ...
                                                    hetomega, hetq, doglobal, fixb, x, p, eps, maxi, betaini, printd, hetdat, hetvar);
        nL = numel(vector_piL);
        Sup_F_piL = NaN(nL, 1); Sup_m_piL = NaN(nL, 1); Sup_pi_piL = NaN(nL, 1); S_piL = cell(nL, 1); nS_piL = NaN(nL, 1);
        for iL = 1:nL
            num_pis = vector_piL(iL):step_pi:1;
            F_ST = zeros(numel(num_pis), max_breaks);
            for m = 1:max_breaks
                for ipi = 1:numel(num_pis)
                    Tpi = floor(num_pis(ipi) * T);
                    if vector_piL(iL) == 1
                        S_m = [1, T];
                    else
                        S_m = S_Dates_m_Heter(y, z, x, datevec, T, num_pis(ipi), m, T_C);
                    end
                    S_index = zeros(T, 1);
                    for irow = 1:size(S_m, 1), S_index(S_m(irow, 1):S_m(irow, 2)) = 1; end
                    y_ST = y(S_index == 1, :); z_ST = z(S_index == 1, :);     % Mx_ST = eye(Tpi) in the legacy code (x empty)
                    F_ST(ipi, m) = F_ST_Stat(y_ST, z_ST, Tpi, q, p, Hac_option);
                end
            end
            [Val_m, ind_m] = max(F_ST, [], 1);
            [Sup_F, Sup_m] = max(Val_m);
            Sup_pi = num_pis(ind_m(Sup_m));
            if vector_piL(iL) == 1, S_sup = [1, T]; else, S_sup = S_Dates_m_Heter(y, z, x, datevec, T, Sup_pi, Sup_m, T_C); end
            Sup_F_piL(iL) = Sup_F; Sup_m_piL(iL) = Sup_m; Sup_pi_piL(iL) = Sup_pi; S_piL{iL} = S_sup;
            nS_piL(iL) = sum(S_sup(:, 2) - S_sup(:, 1) + 1);
            if verbose
                fprintf('  pi_L = %.1f: F* = %.4f at m = %d, pi = %.2f, |S| = %d (Tpi = %d); segments', vector_piL(iL), Sup_F, Sup_m, Sup_pi, nS_piL(iL), floor(Sup_pi*T));
                fprintf(' [%d,%d]', S_sup'); fprintf('\n');
            end
        end
        Results.(D_name) = struct('label', D_list{iD, 2}, 'datevec', datevec, 'Sup_F_piL', Sup_F_piL, 'Sup_m_piL', Sup_m_piL, ...
                                  'Sup_pi_piL', Sup_pi_piL, 'S_piL', {S_piL}, 'nS_piL', nS_piL, 'rej01', Sup_F_piL' > cv01, 'rej05', Sup_F_piL' > cv05);
    end

    txt = fullfile(folder_name, sprintf('Table1_Fstar_h%s_%s.txt', h_rule, timestamp_file));
    fid = fopen(txt, 'w');
    fprintf(fid, 'Table 1: Tests for Identification Failure over all Sub-Populations (F*_T statistic; legacy Run_Fstar_Test_NS18_Heter.m)\n');
    fprintf(fid, 'T=%d, T_C=%d, T_P=%d, m+ = 5, eps1 = 0.05, h = round(eps1*%s) = %d, heteroskedasticity-only variance.\n\n', T, T_C, T_P, h_rule, h);
    fprintf(fid, '%-32s', 'D_t \\ pi_L'); fprintf(fid, ' %12.1f', vector_piL); fprintf(fid, '\n');
    for iD = 1:size(D_list, 1)
        R = Results.(D_list{iD, 1});
        fprintf(fid, '%-32s', D_list{iD, 2}); fprintf(fid, ' %12.2f', R.Sup_F_piL); fprintf(fid, '\n');
        fprintf(fid, '%-32s', '   (m_hat, pi_hat, |S|)'); for iL = 1:numel(vector_piL), fprintf(fid, ' %12s', sprintf('(%d,%.2f,%d)', R.Sup_m_piL(iL), R.Sup_pi_piL(iL), R.nS_piL(iL))); end; fprintf(fid, '\n');
    end
    fprintf(fid, '%-32s', '1%% critical values (Table 4)'); fprintf(fid, ' %12.2f', cv01); fprintf(fid, '\n');
    fprintf(fid, '%-32s', '5%% critical values (Table 4)'); fprintf(fid, ' %12.2f', cv05); fprintf(fid, '\n');
    fprintf(fid, '\nPublished: 30-minute "10^4 x" 95.36 56.50 32.45 18.75 7.42 (scale factor reads 10^3: full-sample statistic 7416); 1-day 155.69 88.22 54.78 33.88 8.09.\n');
    for iD = 1:size(D_list, 1)
        R = Results.(D_list{iD, 1});
        fprintf(fid, '\n%s: S_hat per pi_L\n', D_list{iD, 2});
        for iL = 1:numel(vector_piL), fprintf(fid, '  pi_L = %.1f:', vector_piL(iL)); fprintf(fid, ' [%d,%d]', R.S_piL{iL}'); fprintf(fid, '\n'); end
    end
    fclose(fid);
    save(fullfile(folder_name, sprintf('Table1_Fstar_h%s_%s.mat', h_rule, timestamp_file)), 'Results');
    if verbose, fprintf('\nSaved: %s\n', txt); end
end
