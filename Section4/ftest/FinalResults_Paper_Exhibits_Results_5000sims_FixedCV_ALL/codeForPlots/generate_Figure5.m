% =========================================================================
% GENERATE FIGURE 5 -- Power Curves of F Tests, calibrated NS DGP
% =========================================================================
% Reproduces the aesthetics of the legacy fig1_sec4.m script, reading
% directly from the Paper_Results struct saved by runSection4_Simul.m
% (field naming: Fig5_AR_00, Fig5_AR_25, Fig5_AR_50, Fig5_AR_75).
%
% Expected folder layout (this script lives in .../codeForPlots/):
%   FinalResults_Paper_Exhibits_Results_5000sims_FixedCV_ALL/
%       Paper_Exhibits_Suite_5000sims_ALL_20260704_235754.mat
%       codeForPlots/
%           generate_Figure5.m   <-- this file
% =========================================================================
 
clear; clc; close all;
 
%% 1. Load results
 
mat_path = fullfile('..', 'Paper_Exhibits_Suite_5000sims_ALL_20260704_235754.mat');
load(mat_path, 'Paper_Results');
 
arVec = [0, 0.25, 0.50, 0.75];
n_ar  = length(arVec);
 
thetaVec = [];
FstatTotal           = [];  % rows = arVec, cols = thetaVec -> pi-LATE F*
FstatTotalFullSample  = [];  % rows = arVec, cols = thetaVec -> Full-sample F
 
for k = 1:n_ar
    field = sprintf('Fig5_AR_%02d', arVec(k)*100);
    if ~isfield(Paper_Results, field)
        error('Missing field %s in Paper_Results -- Figure 5 run incomplete.', field);
    end
    res = Paper_Results.(field);
 
    if isempty(thetaVec)
        thetaVec = res.thetaVecPower;
    elseif ~isequal(thetaVec, res.thetaVecPower)
        error(['thetaVecPower differs across arCoef for Figure 5 (field %s). ', ...
               'Cannot stack into a single matrix -- rerun with a consistent grid ', ...
               'across all four arCoef values.'], field);
    end
 
    FstatTotal(k, :)          = res.power_piLATE;      %#ok<SAGROW>
    FstatTotalFullSample(k,:) = res.power_FullSample;  %#ok<SAGROW>
end
 
%% 1b. Restrict theta range for plotting
% The full grid (0:0.05:0.8) includes a region near the DGP pole
% (theta* = T_P/T = 0.075) where D_t is numerically unstable and the
% power curve is non-monotonic (see thetaVecPower <= 0.15). Restricting
% the plotted range to [0.2, 0.75] removes that region and the final
% grid point (0.8) from the figure. This is a plotting-range restriction
% only -- it does not alter the underlying data, and the excluded region
% is still present in Paper_Results / the Paper_Exhibits_*.txt report.
 
theta_plot_min = 0.2;
theta_plot_max = 0.75;
tol = 1e-9;
 
plot_mask = (thetaVec >= theta_plot_min - tol) & (thetaVec <= theta_plot_max + tol);
if ~any(plot_mask)
    error('No grid points fall inside [%.4f, %.4f] -- check thetaVec.', theta_plot_min, theta_plot_max);
end
 
thetaVec_plot              = thetaVec(plot_mask);
FstatTotal_plot             = FstatTotal(:, plot_mask);
FstatTotalFullSample_plot   = FstatTotalFullSample(:, plot_mask);
 
%% 2. Global settings for LaTeX font
 
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
 
%% 3. Generate the combined 2x2 plot
 
fig = figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 10, 8.5]);
sgtitle('Power Curves of $F$ Tests', 'FontWeight', 'bold', 'FontSize', 30, 'Interpreter', 'latex');
 
for k = 1:n_ar
 
    arValue = arVec(k);
 
    ax = subplot(2, 2, k);
    hold(ax, 'on');
 
    powerCurve_Fstar      = FstatTotal_plot(k, :);
    powerCurve_FullSample = FstatTotalFullSample_plot(k, :);
 
    % --- Interpolation for smooth curves (restricted range only) ---
    thetaVecSmooth      = linspace(min(thetaVec_plot), max(thetaVec_plot), 300);
    smooth_Fstar         = interp1(thetaVec_plot, powerCurve_Fstar, thetaVecSmooth, 'spline');
    smooth_FullSample    = interp1(thetaVec_plot, powerCurve_FullSample, thetaVecSmooth, 'spline');
 
    plot(ax, thetaVecSmooth, smooth_Fstar, 'LineStyle', '-', 'Color', [0.0, 0.45, 0.74], 'LineWidth', 3);
    plot(ax, thetaVecSmooth, smooth_FullSample, 'LineStyle', '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 3);
 
    % --- Axes configuration ---
    set(ax, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'XTick', theta_plot_min:0.1:theta_plot_max, ...
        'YTick', 0:0.2:1.0, ...
        'FontSize', 20, ...
        'FontName', 'Computer Modern', ...
        'LineWidth', 0.75);
    set(ax, 'TickLabelInterpreter', 'latex');
 
    axis(ax, [theta_plot_min theta_plot_max 0 1.0]);
 
    title(ax, ['$\rho_e = ' num2str(arValue) '$'], 'FontSize', 24);
    xlabel(ax, '$\theta_{1}$ and $\theta_{3}$', 'FontSize', 20);
 
    if k == 1
        % 'best' = MATLAB's automatic optimal placement (least overlap with
        % plotted data). Legacy used 'SouthEast', which worked for the old
        % [0,0.05] grid where curves rose to the top-left; on the current
        % [0.2,0.75] range the Full-Sample curve occupies that corner.
        % If a fixed position is preferred, 'East' is the closest manual
        % alternative (right side, vertically centered).
        legend(ax, {'$F_{T}^*$ statistic', 'Full sample $F_{T}$ statistic'}, ...
            'Location', 'best', ...
            'FontSize', 14, ...
            'EdgeColor', 'none');
    end
 
    hold(ax, 'off');
end
 
%% 4. Saving (optional)
% output_figure_name = 'Figure5.eps';
% print(fig, output_figure_name, '-depsc');
% fprintf('Figure saved as %s\n', output_figure_name);
 




























































































