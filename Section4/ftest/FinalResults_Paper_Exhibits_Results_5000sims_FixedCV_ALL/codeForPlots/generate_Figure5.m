% =========================================================================
% Generate and Save Combined Power Curve Plot -- FIGURE 5
% =========================================================================
% This script loads simulation results and generates a 2x2 plot comparing
% the power of the F* statistic (subsample-based) and the standard
% full-sample F statistic under different levels of serial correlation.
% The final figure is saved as a high-quality .eps file.
%
% This is a data-source adaptation of the legacy fig1_sec4.m script. The
% ONLY changes relative to the legacy script are:
%   (1) Section 1 ("Setup and Data Loading"): data is now read from the
%       Paper_Results struct saved by runSection4_Simul.m instead of a
%       dedicated 'simulation_results_5000sims_short_sample.mat' file.
%   (2) The plotted theta grid is restricted to the discrete set
%       {0, 0.2, 0.4, 0.6, 0.75} instead of the full simulated grid.
%   (3) XTick is set to [0, 0.4, 0.75] instead of 0:0.025:max(thetaVec).
%   (4) The legend Location is 'best' instead of 'SouthEast'.
% Every other line in Section 2 ("Generate the Combined 2x2 Plot") is
% copied verbatim from fig1_sec4.m: same property order, same FontSize
% values, same explicit per-iteration TickLabelInterpreter reinforcement,
% same reliance (or lack thereof) on groot defaults for title/xlabel.
% =========================================================================
clear;clc;close all;

%% 1. Setup and Data Loading

% --- Load Simulation Results ---
% Paper_Results is the single struct saved by runSection4_Simul.m. It
% contains every scenario (Figure 5, Table 5, Figure 6) as separate
% fields. THIS is the mechanism that guarantees we load Figure 5 data and
% not, e.g., Table 5 or Figure 6 data stored in the same file: the field
% name for each AR(1) coefficient is reconstructed with the exact same
% sprintf pattern used by runSection4_Simul.m to write it
% (sprintf('Fig5_AR_%02d', arCoef*100), e.g. arCoef=0.25 -> 'Fig5_AR_25'),
% and is additionally cross-checked against the internal res.dgp_type
% field ('Sec4'), which runSection4_Simul.m stores inside every Results
% struct it produces. If either check fails, the script stops with an
% explicit error rather than silently plotting the wrong scenario.
mat_path = fullfile('..', 'Paper_Exhibits_Suite_5000sims_ALL_20260704_235754.mat');
load(mat_path, 'Paper_Results');

arVec = [0, 0.25, 0.50, 0.75];

thetaVec_full             = [];
FstatTotal_full           = [];  % rows = arVec, cols = thetaVec_full -> F* power
FstatTotalFullSample_full = [];  % rows = arVec, cols = thetaVec_full -> Full-sample F power

for k = 1:length(arVec)
    field = sprintf('Fig5_AR_%02d', arVec(k)*100);
    if ~isfield(Paper_Results, field)
        error('Missing field %s in Paper_Results -- Figure 5 run incomplete.', field);
    end
    res = Paper_Results.(field);
    if ~strcmp(res.dgp_type, 'Sec4')
        error('Field %s has dgp_type=''%s'', expected ''Sec4''. Data source mismatch.', ...
              field, res.dgp_type);
    end
    if isempty(thetaVec_full)
        thetaVec_full = res.thetaVecPower;
    elseif ~isequal(thetaVec_full, res.thetaVecPower)
        error(['thetaVecPower differs across arCoef for Figure 5 (field %s). ', ...
               'Cannot stack into a single matrix.'], field);
    end
    FstatTotal_full(k, :)           = res.power_piLATE;      %#ok<SAGROW>
    FstatTotalFullSample_full(k, :) = res.power_FullSample;  %#ok<SAGROW>
end

% --- Restrict to the discrete theta grid requested for plotting ---
% Only these five points of the full simulated grid are plotted:
theta_plot_set = [0, 0.2, 0.4, 0.6, 0.75];
tol = 1e-9;

plot_idx = [];
for v = theta_plot_set
    idx = find(abs(thetaVec_full - v) < tol);
    if isempty(idx)
        error('Requested theta value %.4f not found in the simulated grid.', v);
    end
    plot_idx(end+1) = idx; %#ok<SAGROW>
end

thetaVec              = thetaVec_full(plot_idx);
FstatTotal             = FstatTotal_full(:, plot_idx);
FstatTotalFullSample   = FstatTotalFullSample_full(:, plot_idx);

% --- Global Settings for LaTeX Font ---
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');



%% 2. Generate the Combined 2x2 Plot

% --- Create a new figure that will contain all subplots ---
fig = figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 10, 8.5]);

% --- Add the Main Title for the Entire Figure ---
sgtitle('Power Curves of $F$ Tests', 'FontWeight', 'bold', 'FontSize', 30, 'Interpreter','Latex');

% --- Loop through each AR coefficient to create a separate subplot ---
for k = 1:length(arVec)

    arValue = arVec(k);

    % --- Create a new subplot in a 2x2 grid ---
    ax = subplot(2, 2, k);
    hold on;

    % --- Extract the data for the current AR coefficient ---
    powerCurve_Fstar = FstatTotal(k, :);
    powerCurve_FullSample = FstatTotalFullSample(k, :);

    % --- Interpolation for smooth curves ---
    % NOTE: pchip (not spline) -- with only 5 sparse grid points and a
    % steep transition between theta=0.2 and theta=0.4, a natural cubic
    % spline overshoots well outside [0,1] (verified: down to -0.21, up
    % to 1.05), which gets clipped by the axis limits and renders as a
    % broken/garbled curve. pchip is shape-preserving and guarantees the
    % interpolant never exceeds the local min/max of the data points.
    thetaVecSmooth = linspace(min(thetaVec), max(thetaVec), 300);
    smooth_Fstar = interp1(thetaVec, powerCurve_Fstar, thetaVecSmooth, 'pchip');
    smooth_FullSample = interp1(thetaVec, powerCurve_FullSample, thetaVecSmooth, 'pchip');

    % --- Plot the two power curves with specified styles ---
    plot(thetaVecSmooth, smooth_Fstar, 'LineStyle', '-', 'Color', [0.0, 0.45, 0.74], 'LineWidth', 3);
    plot(thetaVecSmooth, smooth_FullSample, 'LineStyle', '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 3);

    % --- Axes Configuration (No Grid) ---
    set(ax, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'XTick', [0, 0.4, 0.75], ...
        'YTick', 0:0.2:1.0, ...
        'FontSize', 20, ...
        'FontName', 'Computer Modern', ...
        'LineWidth', 0.75);
    set(gca,'TickLabelInterpreter','latex')
    % --- Set axis limits ---
    axis([min(thetaVec) max(thetaVec) 0 1.0]);

    % --- Subplot Title ---
    title(['$\rho_e = ' num2str(arValue) '$'], 'FontSize', 24);

    % --- Add Axis Labels ---
    % Add x-label to all plots, as in the original
    xlabel('$\theta_{1}$ and $\theta_{3}$', 'FontSize', 20);
    % --- Add a legend to the first subplot only for clarity ---
    if k == 1
        legend({'$F_{T}^*$ statistic', 'Full sample $F_{T}$ statistic'}, ...
            'Location', 'best', ...
            'FontSize', 14, ...
            'EdgeColor', 'none');
    end
    hold off;
end

% --- 4. Saving (optional) ---
% output_figure_name = 'Figure5.eps';
% print(fig, output_figure_name, '-depsc');
% fprintf('Figure saved as %s\n', output_figure_name);