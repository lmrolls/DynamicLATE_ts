% =========================================================================
% MATLAB SCRIPT TO GENERATE A 2x2 SUBPLOT OF F-TEST POWER CURVES -- FIGURE 6
% =========================================================================
% This script creates a single figure with 4 subplots in a 2x2 grid.
% It adopts the formatting from 'fig1_sec4.m' and combines it with a
% row/column labeling strategy seen in other figures.
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
% This is a data-source adaptation of the legacy Figure_2_Sect_4.m script.
% The ONLY change relative to the legacy script is in Section 3
% ("Looping and Plotting"): instead of loading a dedicated .mat file per
% (pi0, rho) combination via isfile/load, each combination is now read as
% a field of the single Paper_Results struct saved by runSection4_Simul.m.
% The field name is reconstructed with the exact same sprintf pattern
% runSection4_Simul.m uses to write it (sprintf('Fig6_pi0_%02d_AR_%02d',
% pi0*10, rho*100), e.g. pi0=0.6, rho=0.25 -> 'Fig6_pi0_06_AR_25'), and is
% additionally cross-checked against the internal res.dgp_type field
% ('M1_reverse') and the stored res.pi0/res.arCoef, which
% runSection4_Simul.m writes inside every Results struct it produces.
% The graceful skip-with-placeholder behavior of the legacy script (for a
% missing file) is preserved, now keyed on a missing struct field instead
% of a missing file. No other line is changed: same property order, same
% FontSize values, same title/xlabel/legend arguments, same row-label
% mechanism via hidden axes.
% =========================================================================

clear; clc; close all;

% --- 1. CONFIGURATION ---

% --- Global Settings for LaTeX Font (from fig1_sec4.m) ---
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% --- Load the unified results struct produced by runSection4_Simul.m ---
% Paper_Results contains every scenario (Figure 5, Table 5, Figure 6) as
% separate fields; only the Fig6_pi0_XX_AR_XX fields are used below.
mat_path = fullfile('..', 'Paper_Exhibits_Suite_5000sims_ALL_20260704_235754.mat');
load(mat_path, 'Paper_Results');

% --- Define the parameter values for the labels ---
rho_values = [0.25, 0.75];
pi0_values = [0.6, 0.8];

% --- T is not stored per-scenario in Results for the M1/M2 branch; it is
% hardcoded to 200 inside Run_MC_Scenario (runSection4_Simul.m) and never
% propagated to the saved struct. Hardcoded here to match. ---
T_M1 = 200;

% --- 2. FIGURE CREATION ---

fig = figure('Color', 'white', 'Units', 'inches', 'Position', [1, 1, 10, 9.5]);
% --- Add the Main Title for the Entire Figure (using sgtitle from fig1_sec4.m) ---
sgtitle('Power Curves of $F$ Tests', 'Interpreter', 'latex', 'FontSize', 30);

% --- 3. LOOPING AND PLOTTING ---

for row = 1:2
    for col = 1:2

        ax = subplot(2, 2, (row-1)*2 + col); % Use subplot directly for better control with sgtitle
        all_axes(row, col) = ax; % Store handle
        hold(ax, 'on');

        pi0 = pi0_values(row);
        rho = rho_values(col);
        field = sprintf('Fig6_pi0_%02d_AR_%02d', pi0*10, rho*100);

        if ~isfield(Paper_Results, field)
            fprintf('Warning: Field not found, skipping subplot. -> %s\n', field);
            text(0.5, 0.5, 'Data Not Found', 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
            axis off;
            continue; % Skip to the next iteration
        end

        fprintf('Processing: %s\n', field);
        res = Paper_Results.(field);
        if ~strcmp(res.dgp_type, 'M1_reverse') || res.pi0 ~= pi0 || res.arCoef ~= rho
            error('Field %s does not match expected (dgp_type=M1_reverse, pi0=%.2f, arCoef=%.2f). Data source mismatch.', ...
                  field, pi0, rho);
        end

        % Load data, including dVecPower and T for the new x-axis
        dVecPower                = res.thetaVecPower;
        T                        = T_M1;
        FstatTotal_SA_piLATE     = res.power_piLATE;
        FstatTotal_SA_FullSample = res.power_FullSample;

        % --- Calculate the new horizontal axis values ---
        x_axis_values = dVecPower / sqrt(T);

        % --- Interpolation for smooth curves using the new x-axis ---
        x_axis_smooth = linspace(min(x_axis_values), max(x_axis_values), 300);
        smooth_Fstar = interp1(x_axis_values, FstatTotal_SA_piLATE, x_axis_smooth, 'spline');
        smooth_FullSample = interp1(x_axis_values, FstatTotal_SA_FullSample, x_axis_smooth, 'spline');

        % --- Plot the two power curves with specified styles ---
        % F* (piLATE) is solid blue
        plot(ax, x_axis_smooth, smooth_Fstar, 'LineStyle', '-', 'Color', [0.0, 0.45, 0.74], 'LineWidth', 3);
        % F (Full Sample) is dashed red
        plot(ax, x_axis_smooth, smooth_FullSample, 'LineStyle', '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 3);

        % --- Axes Configuration (from fig1_sec4.m) ---
        set(ax, ...
            'Box', 'on', ...
            'TickDir', 'in', ...
            'XTick', 0:0.5:ceil(max(x_axis_values)), ... % Adjusted for new axis
            'YTick', 0:0.2:1.0, ...
            'FontSize', 18, ...
            'FontName', 'Computer Modern', ...
            'LineWidth', 0.75);

        % --- Set axis limits for the new x-axis ---
        axis(ax, [min(x_axis_values) max(x_axis_values) 0 1.0]);

        % --- Add a title to every subplot (format from fig1_sec4.m) ---
        title(ax, ['$\rho_e = \rho_u = ' num2str(rho_values(col)) '$'], 'Interpreter', 'latex', 'FontSize', 20);

        % --- Add new Axis Label ---
        xlabel(ax, '$d/\sqrt{T}$', 'FontSize', 24);

        % --- Add a legend to the first subplot only (from fig1_sec4.m) ---
        if row == 1 && col == 1
            legend(ax, {'$F_{T}^*$ statistic', 'Full sample $F_{T}$ statistic'}, ...
                'Location', 'SouthEast', ...
                'FontSize', 12, ...
                'EdgeColor', 'none');
        end

        hold(ax, 'off');
    end
end

% --- 4. ROW LABELS ---
% Add the 'pi_0=...' labels BELOW each row using a hidden axes.
ax_hidden = axes('Position',[0 0 1 1],'Visible','off');
drawnow; % Ensure layout is rendered before getting positions

for row = 1:2
    pos_left = get(all_axes(row, 1), 'Position');
    pos_right = get(all_axes(row, 2), 'Position');

    % Calculate the horizontal center between the two plots
    x_center = pos_left(1) + 0.5 * (pos_right(1) + pos_right(3) - pos_left(1)) - 0.008;

    % Calculate the vertical position BELOW the plots
    % NOTE: This offset may need minor adjustment for perfect alignment.
    y_pos = pos_left(2) + 0.365;

    text(ax_hidden, x_center, y_pos, sprintf('$\\pi_0 = %.1f$', pi0_values(row)), ...
               'Interpreter', 'latex', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'top', ...
               'FontSize', 24, 'FontWeight', 'normal');
end

% --- 5. SAVING (OPTIONAL) ---
% output_figure_name = 'figure_section4_new.png';
% saveas(fig, output_figure_name);
% fprintf('Figure saved as %s\n', output_figure_name);