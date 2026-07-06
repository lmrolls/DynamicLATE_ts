% =========================================================================
% GENERATE FIGURE 6 -- Power Curves of F Tests, proportional anomaly (M1)
% =========================================================================
% Reproduces the aesthetics of the legacy Figure_2_Sect_4.m script,
% reading directly from the Paper_Results struct saved by
% runSection4_Simul.m (field naming: Fig6_pi0_06_AR_25, etc.).
%
% LAYOUT:
%                  rho_e=rho_u = 0.25    rho_e=rho_u = 0.75
%               +---------------------+---------------------+
%               | Legend Here         |                     |
%               | (1,1)               | (1,2)               |
%               +---------------------+---------------------+
%                             pi_0 = 0.6
%               +---------------------+---------------------+
%               |                     |                     |
%               | (2,1)               | (2,2)               |
%               +---------------------+---------------------+
%                             pi_0 = 0.8
%
% Expected folder layout (this script lives in .../codeForPlots/):
%   FinalResults_Paper_Exhibits_Results_5000sims_FixedCV_ALL/
%       Paper_Exhibits_Suite_5000sims_ALL_20260704_235754.mat
%       codeForPlots/
%           generate_Figure6.m   <-- this file
% =========================================================================

clear; clc; close all;

%% 1. Load results

mat_path = fullfile('..', 'Paper_Exhibits_Suite_5000sims_ALL_20260704_235754.mat');
load(mat_path, 'Paper_Results');

% NOTE: T is not stored in Results for the M1/M2 branch in
% runSection4_Simul.m -- it is hardcoded to 200 inside Run_MC_Scenario and
% never propagated to the saved struct. Hardcoded here to match. If T is
% ever changed for M1/M2, update this value (or, better, add
% Results.T = T; inside Run_MC_Scenario so this stops being manual).
T_M1 = 200;

pi0_values  = [0.6, 0.8];
rho_values  = [0.25, 0.75];

%% 2. Global settings for LaTeX font

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% 3. Figure creation

fig = figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 10, 9.5]);
sgtitle('Power Curves of $F$ Tests', 'Interpreter', 'latex', 'FontSize', 30);

all_axes = gobjects(2, 2);
first_dvec_grid = [];

%% 4. Looping and plotting

for row = 1:2
    for col = 1:2

        pi0 = pi0_values(row);
        rho = rho_values(col);

        ax = subplot(2, 2, (row-1)*2 + col);
        all_axes(row, col) = ax;
        hold(ax, 'on');

        field = sprintf('Fig6_pi0_%02d_AR_%02d', pi0*10, rho*100);
        if ~isfield(Paper_Results, field)
            error('Missing field %s in Paper_Results -- Figure 6 run incomplete.', field);
        end
        res = Paper_Results.(field);

        dVecPower                = res.thetaVecPower;    % field kept as-is in new struct; semantically "d" here
        FstatTotal_SA_piLATE     = res.power_piLATE;
        FstatTotal_SA_FullSample = res.power_FullSample;

        if isempty(first_dvec_grid)
            first_dvec_grid = dVecPower;
        elseif ~isequal(first_dvec_grid, dVecPower)
            warning(['dVecPower differs across (pi0, rho) combinations for Figure 6 ', ...
                     '(field %s). Each subplot uses its own x-axis independently, so ', ...
                     'this does not break the plot, but check whether this is intended.'], field);
        end

        % --- Horizontal axis: d / sqrt(T) ---
        x_axis_values = dVecPower / sqrt(T_M1);

        % --- Interpolation for smooth curves ---
        x_axis_smooth      = linspace(min(x_axis_values), max(x_axis_values), 300);
        smooth_Fstar        = interp1(x_axis_values, FstatTotal_SA_piLATE, x_axis_smooth, 'spline');
        smooth_FullSample   = interp1(x_axis_values, FstatTotal_SA_FullSample, x_axis_smooth, 'spline');

        plot(ax, x_axis_smooth, smooth_Fstar, 'LineStyle', '-', 'Color', [0.0, 0.45, 0.74], 'LineWidth', 3);
        plot(ax, x_axis_smooth, smooth_FullSample, 'LineStyle', '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 3);

        set(ax, ...
            'Box', 'on', ...
            'TickDir', 'in', ...
            'XTick', 0:0.5:ceil(max(x_axis_values)), ...
            'YTick', 0:0.2:1.0, ...
            'FontSize', 18, ...
            'FontName', 'Computer Modern', ...
            'LineWidth', 0.75);

        axis(ax, [min(x_axis_values) max(x_axis_values) 0 1.0]);

        title(ax, ['$\rho_e = \rho_u = ' num2str(rho) '$'], 'Interpreter', 'latex', 'FontSize', 20);
        xlabel(ax, '$d/\sqrt{T}$', 'FontSize', 24);

        if row == 1 && col == 1
            legend(ax, {'$F_{T}^*$ statistic', 'Full sample $F_{T}$ statistic'}, ...
                'Location', 'SouthEast', ...
                'FontSize', 12, ...
                'EdgeColor', 'none');
        end

        hold(ax, 'off');
    end
end

%% 5. Row labels: pi_0 = ... below each row

ax_hidden = axes('Position', [0 0 1 1], 'Visible', 'off');
drawnow;

for row = 1:2
    pos_left  = get(all_axes(row, 1), 'Position');
    pos_right = get(all_axes(row, 2), 'Position');

    x_center = pos_left(1) + 0.5 * (pos_right(1) + pos_right(3) - pos_left(1)) - 0.008;
    y_pos    = pos_left(2) + 0.365;

    text(ax_hidden, x_center, y_pos, sprintf('$\\pi_0 = %.1f$', pi0_values(row)), ...
        'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 24, 'FontWeight', 'normal');
end

%% 6. Saving (optional)
% output_figure_name = 'Figure6.png';
% saveas(fig, output_figure_name);
% fprintf('Figure saved as %s\n', output_figure_name);
