function Results = runSection7_Figure2(external_dir, break_dates, verbose)
% RUNSECTION7_FIGURE2  Figure 2 of the paper (Section 7.1): 1-day changes in 2-Year Treasury yields in the control
% sample (top) and in the policy sample (bottom), with the break dates estimated in the control sample by the
% Casini-Perron (2024) "All-Inside" change-point procedure for the spectrum.
%
% Legacy: pi-LATE\Figure_Breaks_CasiniPerron2024\All-Inside_Matlab\All_Inside.m (Casini-Perron code, 116 files;
% runs on DNY_C from nak_data.mat, saves Break_Dates_Final in Breakdates.mat; published dates [320 436 583] in the
% control sample of 762 days) and Figure1.m of that folder (plot). The All-Inside code is NOT part of this
% repository: it is an external dependency, like Lewis's code for Table 2. Copy the folder All-Inside_Matlab into
% external_dir (default Section7\external\AllInside, kept out of git) or pass its path:
%   * if All_Inside.m is found there, the procedure is run (cwd = that folder, as the script expects its data
%     nak_data.mat and sub-folders) and the break dates are read from the Breakdates.mat it saves;
%   * else, if Breakdates.mat is found there, its Break_Dates_Final is used;
%   * else the break dates given in break_dates (default: the published [320 436 583]) are plotted and flagged.
% The series plotted are DNY2 (observations 401:end of ..\m_break\dataForMatlabAltEstTick.csv) split by FOMCused,
% identical to DNY_C / DNY_P of nak_data.mat. Output: .txt with the break dates and their calendar dates, .png.
% =========================================================================
    if nargin < 1 || isempty(external_dir), external_dir = fullfile('external', 'AllInside'); end
    if nargin < 2 || isempty(break_dates), break_dates = [320 436 583]; end
    if nargin < 3 || isempty(verbose), verbose = true; end
    folder_name = 'Paper_Exhibits_Results_SEC7_FIGURE2';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    timestamp_file = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    here = pwd;

    dta = readtable(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv'));
    dates = read_first_column(fullfile('..', 'm_break', 'dataForMatlabAltEstTick.csv')); dates = dates(401:end);
    FOMC = dta.FOMCused(401:end); DNY2 = dta.DNY2(401:end);
    DNY_C = DNY2(FOMC == 0); DNY_P = DNY2(FOMC == 1); dates_C = dates(FOMC == 0); dates_P = dates(FOMC == 1);

    source = 'published values (All-Inside code not found)';
    if exist(fullfile(external_dir, 'All_Inside.m'), 'file') == 2
        try
            cd(external_dir); evalc('run(''All_Inside.m'')'); cd(here);
            L = load(fullfile(external_dir, 'Breakdates.mat')); break_dates = double(L.Break_Dates_Final(:)'); source = 'All_Inside.m run now';
        catch ME
            cd(here); if verbose, fprintf('All_Inside.m failed: %s\n', ME.message); end
            if exist(fullfile(external_dir, 'Breakdates.mat'), 'file') == 2
                L = load(fullfile(external_dir, 'Breakdates.mat')); break_dates = double(L.Break_Dates_Final(:)'); source = 'Breakdates.mat in external_dir';
            end
        end
    elseif exist(fullfile(external_dir, 'Breakdates.mat'), 'file') == 2
        L = load(fullfile(external_dir, 'Breakdates.mat')); break_dates = double(L.Break_Dates_Final(:)'); source = 'Breakdates.mat in external_dir';
    end
    Results = struct('break_dates', break_dates, 'break_calendar', {dates_C(break_dates)}, 'source', source, 'T_C', numel(DNY_C), 'T_P', numel(DNY_P));
    txt = fullfile(folder_name, sprintf('Figure2_Breaks_%s.txt', timestamp_file)); fid = fopen(txt, 'w');
    fprintf(fid, 'Figure 2: break dates in the control sample of 2-Year yields (Casini-Perron 2024 All-Inside). Source: %s\n', source);
    fprintf(fid, 'Indices (control sample, T_C = %d): %s; calendar dates: %s\n', numel(DNY_C), mat2str(break_dates), strjoin(dates_C(break_dates)', ', '));
    fprintf(fid, 'Published: [320 436 583].\n'); fclose(fid);
    if verbose, fprintf('Break dates (%s): %s = %s\n', source, mat2str(break_dates), strjoin(dates_C(break_dates)', ', ')); end
    try
        fh = figure('Visible', 'off');
        subplot(2, 1, 1); plot(1:numel(DNY_C), DNY_C, 'LineWidth', 1.2); hold on;
        for b = break_dates, xline(b, '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.5); end
        xlim([1 numel(DNY_C)]); title('2-Year Treasury yields (1-day changes) and structural breaks: control sample');
        tk = round(linspace(1, numel(DNY_C), 6)); set(gca, 'XTick', tk, 'XTickLabel', dates_C(tk));
        subplot(2, 1, 2); plot(1:numel(DNY_P), DNY_P, 'LineWidth', 1.2); xlim([1 numel(DNY_P)]); title('Policy sample');
        tk = round(linspace(1, numel(DNY_P), 6)); set(gca, 'XTick', tk, 'XTickLabel', dates_P(tk));
        print(fh, fullfile(folder_name, sprintf('Figure2_Breaks_%s.png', timestamp_file)), '-dpng'); close(fh);
    catch ME
        if verbose, fprintf('Figure 2 not drawn: %s\n', ME.message); end
    end
    save(fullfile(folder_name, sprintf('Figure2_Breaks_%s.mat', timestamp_file)), 'Results');
end

function c = read_first_column(fname)
    fid = fopen(fname, 'r'); fgetl(fid); c = {};
    while true
        l = fgetl(fid); if ~ischar(l), break; end
        k = find(l == ',', 1); if isempty(k), c{end+1, 1} = strtrim(l); else, c{end+1, 1} = strtrim(l(1:k-1)); end
    end
    fclose(fid);
end
