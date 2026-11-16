function Results = runSection7_Figure1(verbose)
% RUNSECTION7_FIGURE1  Figure 1 of the paper (Section 1): distributions of the full-sample first-stage F and of the
% F* statistics in the survey of IV applications (AER, Econometrica, JPE, QJE, REStud, 2019-2022), time series
% specifications (top panels) and panel data specifications (bottom panels); and the summary numbers quoted in the
% text (1,560 specifications from 18 papers: 199 time series, 1,361 panel; 75% / 72% with F below 24; median F 12.63 / 9.29).
%
% Legacy script: pi-LATE\Figure_Historgram\Figure_Histogram_time_series_panel.m (7/7/2026) with the survey data
% Data_Survey_top5_time_series.xlsx and Data_Survey_top5_Panel.xlsx (15/9/2025; copies in Section7\data\; columns:
% full-sample F, weak-identification dummy, F*, pi_hat). As in the script: histograms with bin edges 0:5:100 after
% dropping the specifications with F* > 100 and then those with F > 100; vertical lines at 3.841 and 6.34 for F
% (the script's values; chi2(1) critical values are 3.841 and 6.635) and at 12.74 and 17.89 for F* (Table 4, q = 1,
% pi_L = 0.6, 5% and 1%; the published caption says "8.28 and 11.63"). The script's statistics (mean, median,
% quartiles, share of specifications flagged as weakly identified, share of those with F* >= 24) are reported.
% The figure is saved as .png (plotting guarded).
% =========================================================================
    if nargin < 1 || isempty(verbose), verbose = true; end
    folder_name = 'Paper_Exhibits_Results_SEC7_FIGURE1';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    files = {'Data_Survey_top5_time_series.xlsx', 'Time series'; 'Data_Survey_top5_Panel.xlsx', 'Panel data'};
    Results = struct();
    txt = fullfile(folder_name, sprintf('Figure1_Survey_%s.txt', timestamp_file)); fid = fopen(txt, 'w');
    fprintf(fid, 'Figure 1 / survey summary (legacy Figure_Histogram_time_series_panel.m; data Section7\\data\\*.xlsx)\n');
    for k = 1:2
        M = read_survey(fullfile('data', files{k, 1}));
        F = M(:, 1); weak = M(:, 2); Fstar = M(:, 3); pih = M(:, 4);
        S = struct('n', numel(F), 'mean_F', mean(F), 'median_F', median(F), 'q25_F', quantile(F, 0.25), 'q75_F', quantile(F, 0.75), ...
                   'share_F_below_24', mean(F < 24), 'share_F_below_10', mean(F < 10), 'share_weak_dummy', mean(weak), ...
                   'mean_Fstar', mean(Fstar), 'median_Fstar', median(Fstar), 'q25_Fstar', quantile(Fstar, 0.25), 'q75_Fstar', quantile(Fstar, 0.75), ...
                   'percent_strong_among_weak', 100 * mean(Fstar(weak == 1) >= 24), 'mean_pihat', mean(pih), 'median_pihat', median(pih), ...
                   'q25_pihat', quantile(pih, 0.25), 'q75_pihat', quantile(pih, 0.75));
        keep = Fstar <= 100; keep = keep & F <= 100;      % script: drop F* > 100, then F > 100
        S.F_cut = F(keep); S.Fstar_cut = Fstar(keep); S.n_cut = nnz(keep);
        Results.(strrep(files{k, 2}, ' ', '_')) = S;
        line = sprintf('%s: n = %d; F: mean %.2f median %.2f [q25 %.2f, q75 %.2f], share < 24 = %.3f, share < 10 = %.3f, weak dummy %.3f; F*: mean %.2f median %.2f [%.2f, %.2f], F* >= 24 among weak %.1f%%; pi_hat: mean %.3f median %.2f [%.2f, %.2f]; plotted (F*, F <= 100): %d', ...
            files{k, 2}, S.n, S.mean_F, S.median_F, S.q25_F, S.q75_F, S.share_F_below_24, S.share_F_below_10, S.share_weak_dummy, S.mean_Fstar, S.median_Fstar, S.q25_Fstar, S.q75_Fstar, S.percent_strong_among_weak, S.mean_pihat, S.median_pihat, S.q25_pihat, S.q75_pihat, S.n_cut);
        if verbose, fprintf('%s\n', line); end
        fprintf(fid, '%s\n', line);
    end
    fprintf(fid, 'Published text: 1,560 specifications (199 time series, 1,361 panel); 75%% / 72%% with F below 24; median F 12.63 / 9.29.\n');
    fclose(fid);
    save(fullfile(folder_name, sprintf('Figure1_Survey_%s.mat', timestamp_file)), 'Results');
    try
        fh = figure('Visible', 'off'); edges = 0:5:100; k = 0;
        for nm = {'Time_series', 'Panel_data'}
            S = Results.(nm{1}); k = k + 1;
            subplot(2, 2, 2*k - 1); histogram(S.F_cut, 'BinEdges', edges, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none'); hold on;
            xline(3.841, '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 2); xline(6.34, '-', 'Color', [0.735 0.078 0.084], 'LineWidth', 2);
            title(sprintf('%s: full sample first-stage F', strrep(nm{1}, '_', ' ')));
            subplot(2, 2, 2*k); histogram(S.Fstar_cut, 'BinEdges', edges, 'FaceColor', [0.466 0.674 0.388], 'EdgeColor', 'none'); hold on;
            xline(12.74, '-', 'Color', [0.85 0.325 0.098], 'LineWidth', 2); xline(17.89, '-', 'Color', [0.735 0.078 0.084], 'LineWidth', 2);
            title(sprintf('%s: F* statistic', strrep(nm{1}, '_', ' ')));
        end
        print(fh, fullfile(folder_name, sprintf('Figure1_Survey_%s.png', timestamp_file)), '-dpng'); close(fh);
    catch ME
        if verbose, fprintf('Figure 1 not drawn: %s\n', ME.message); end
    end
    if verbose, fprintf('Saved: %s\n', txt); end
end

function M = read_survey(f)
    if exist('readmatrix', 'file') == 2
        M = readmatrix(f);
    else
        M = xlsread(f);
    end
    M = M(all(~isnan(M(:, 1:4)), 2), 1:4);
end
