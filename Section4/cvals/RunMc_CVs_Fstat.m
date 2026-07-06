function RunMc_CVs_Fstat()
% RUNMC_CVS_FSTAT Suite automatizzata per generare le tabelle dei CV asintotici.
% Salva i risultati in formato .mat e in un file .txt formattato LaTeX.

    clear; clc;
    restoredefaultpath;
    addpath('..\..\m_break\');

    %% 1. Configurazione della Griglia
    q_values      = [1, 2, 3, 4, 5, 10];
    piLow_values  = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    T             = 1000;
    totreps       = 10000; 

    num_pi = length(piLow_values);
    num_q  = length(q_values);

    Table_CV_90 = zeros(num_pi, num_q);
    Table_CV_95 = zeros(num_pi, num_q);
    Table_CV_99 = zeros(num_pi, num_q);

    row_names = arrayfun(@(x) sprintf('%.2f', x), piLow_values, 'UniformOutput', false);
    col_names = arrayfun(@(x) sprintf('%d', x), q_values, 'UniformOutput', false);

    fprintf('============================================================\n');
    fprintf('AVVIO REPLICAZIONE TABELLE CV ASINTOTICI\n');
    fprintf('============================================================\n\n');

    %% 2. Loop Principale
    for i = 1:num_pi
        for j = 1:num_q
            fprintf('Calcolo pi_L=%.2f, q=%d... ', piLow_values(i), q_values(j));
            results = Compute_CVs(q_values(j), piLow_values(i), T, totreps, [], true);
            
            Table_CV_90(i, j) = results.CV_90;
            Table_CV_95(i, j) = results.CV_95;
            Table_CV_99(i, j) = results.CV_99;
            fprintf('Fatto.\n');
        end
    end

    %% 3. Gestione Output e Salvataggio
    folder_name = 'Results_CV_Export';
    if ~exist(folder_name, 'dir'), mkdir(folder_name); end
    
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    mat_filename = fullfile(folder_name, sprintf('CV_Results_%s.mat', timestamp));
    save(mat_filename, 'Table_CV_90', 'Table_CV_95', 'Table_CV_99', 'q_values', 'piLow_values');
    
    % Generazione file TXT con output LaTeX
    txt_filename = fullfile(folder_name, sprintf('CV_Tables_LaTeX_%s.txt', timestamp));
    fid = fopen(txt_filename, 'w');
    
    write_latex_table(fid, '10% Critical Values', Table_CV_90, row_names, col_names);
    write_latex_table(fid, '5% Critical Values',  Table_CV_95, row_names, col_names);
    write_latex_table(fid, '1% Critical Values',  Table_CV_99, row_names, col_names);
    
    fclose(fid);
    fprintf('\nEsecuzione completata. Dati e LaTeX salvati in: %s\n', folder_name);
end

%% ------------------------------------------------------------------------
%  Helper Function: Scrittura formato LaTeX
%  ------------------------------------------------------------------------
function write_latex_table(fid, title_str, data, rows, cols)
    fprintf(fid, '%% %s\n', title_str);
    fprintf(fid, '\\begin{table}[h]\n\\centering\n\\begin{tabular}{l|%s}\n', repmat('c', 1, length(cols)));
    fprintf(fid, '$\\pi_L \\setminus q$ & %s \\\\ \\hline\n', strjoin(cols, ' & '));
    
    for r = 1:length(rows)
        row_str = sprintf('%.2f & ', str2double(rows{r}));
        row_data = sprintf('%.2f & ', data(r, :));
        fprintf(fid, '%s %s \\\\\n', row_str, row_data(1:end-2));
    end
    
    fprintf(fid, '\\end{tabular}\n\\caption{%s}\n\\end{table}\n\n', title_str);
end