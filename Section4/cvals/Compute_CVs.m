function [CV_Table] = Compute_CVs(q_value, piLow, T, totreps, qw, set_seed)
% COMPUTE_CVS Calcola i Valori Critici per il test Dynamic LATE
%
% INPUTS:
%   q_value   : (Obbligatorio) Dimensione del processo (numero di strumenti)
%   piLow     : (Opzionale) Limite inferiore per pi (Default: 0.6)
%   T         : (Opzionale) Dimensione campionaria temporale (Default: 1000)
%   totreps   : (Opzionale) Numero di replicazioni Monte Carlo (Default: 10000)
%   qw        : (Opzionale) Sovrascrive q_value se diverso (Default: = q_value)
%   set_seed  : (Opzionale) Booleano, se true imposta RNG deterministico (Default: true)
%
% OUTPUT:
%   CV_Table  : Tabella con i valori critici al 90%, 95%, 99%

    %% 1. Gestione Input e Defaults
    if nargin < 1
        error('Il parametro q_value è obbligatorio.');
    end

    if nargin < 2 || isempty(piLow)
        piLow = 0.6;
    end

    if nargin < 3 || isempty(T)
        T = 1000;
    end

    if nargin < 4 || isempty(totreps)
        totreps = 10000;
    end

    if nargin < 5 || isempty(qw)
        qw = q_value; % Di default qw è uguale al q_value passato
    end

    if nargin < 6 || isempty(set_seed)
        set_seed = true;
    end

    %% 2. Impostazioni Parametriche (Hardcoded come da richiesta)
    % Impostazioni per Procedure_piL_Full_Sample
    Z = ones(T-1,1);
    X = [];
    q_procedure = size(Z,2); % Nota: questo è diverso da qw
    p = size(X,2);

    m_breaks = 5; 
    eps1 = .05; 
    h = round(eps1*T);    
    fixb=0; betaini=0; maxi=20; printd=0; eps=0.0001; robust=1; 
    prewhit=0; hetdat=1; hetvar=1; hetomega=1; hetq=1; doglobal=1;

    % Impostazioni griglia Pi
    vector_piL = piLow:0.1:1; 
    step_pi = 0.05;
    max_breaks = 5;

    fprintf('Avvio simulazione CVs: q=%d, piLow=%.2f, T=%d, Reps=%d\n', qw, piLow, T, totreps);

    %% 3. Simulazione Monte Carlo
    outReps = zeros(totreps, 1);

    % Costruzione dell'offset algebrico per la separazione degli stream stocastici.
    % La spaziatura garantisce una distanza minima di 100.000 tra seed di griglie adiacenti,
    % prevenendo la sovrapposizione delle traiettorie per totreps = 10000.
    offset = q_value * 1e7 + round(piLow * 100) * 1e5 + T * 1e2 + totreps; 

    parfor it = 1:totreps
        % Gestione Seed indicizzata per iterazione e configurazione parametrica
        if set_seed
            rng(it + offset, 'twister');
        end

        % Generazione Processo di Wiener Multivariato vettorializzato in O(T)
        D = cumsum(randn(T, qw)) ./ sqrt(T);

        % Isolamento del primo operatore di differenza per il feed in Bai-Perron
        fd_D = D(2:end, :) - D(1:end-1, :);

        % Inizializzazione contenitore Supremo
        Sup_F_qw = zeros(qw, 1);

        % Loop sui k strumenti (dimensioni del Wiener)
        for k = 1:qw
            % NOTA: Procedure_piL_Full_Sample deve essere nel path di MATLAB
            [~,datevec,~] = Procedure_piL_Full_Sample(T-1, fd_D(:,k), Z, q_procedure, m_breaks, h, eps1, robust, prewhit, hetomega, hetq, doglobal, fixb, X, p, eps, maxi, betaini, printd, hetdat, hetvar);

            % Correzione indici datevec
            for ii = 1:size(datevec,2)
                datevec(1:ii,ii) = datevec(1:ii,ii) + 1;
            end

            Sup_F_piL = zeros(length(vector_piL), 1);

            % Loop su vari pi_Lower nella griglia (se vector_piL ha più elementi)
            for piL_idx = 1:length(vector_piL)
                num_pis = vector_piL(piL_idx):step_pi:1;
                index_pis = length(num_pis);
                F_S = zeros(index_pis, max_breaks);

                for m = 1:max_breaks
                    for pi_s_idx = 1:index_pis
                        current_pi = num_pis(pi_s_idx);

                       
                        S_m = S_Dates_Norm(D, datevec, T, current_pi, m);

                        if isempty(S_m)
                            error('Error: S_m is empty');
                        end

                        S_m = sortrows(S_m, 1); 

                        if size(S_m,1) > 1
                            [S_m] = merge_Dates(S_m); 
                        end

                        vec_sum = zeros(1, qw);
                        for i = 1:size(S_m,1)
                            if S_m(i,1) == 1
                                vec_sum = vec_sum + D(S_m(i,2),:);
                            else
                                vec_sum = vec_sum + (D(S_m(i,2),:) - D(S_m(i,1) - 1,:));
                            end
                        end
                        out = norm(vec_sum)^2 / (qw * current_pi);
                        F_S(pi_s_idx, m) = out;
                    end
                end
                [Sup_F, ~] = max(F_S, [], 'all');
                Sup_F_piL(piL_idx) = Sup_F; 
            end

            [Sup_F, ~] = max(Sup_F_piL, [], 'all');
            Sup_F_qw(k,1) = Sup_F;
        end

        [value, ~] = max(Sup_F_qw);
        outReps(it) = value;
    end

    %% 4. Calcolo Quantili e Output
    ss = sort(outReps);
    cv_90 = ss(round(totreps * 0.90));
    cv_95 = ss(round(totreps * 0.95));
    cv_99 = ss(round(totreps * 0.99));

    CV_Table = table(cv_90, cv_95, cv_99, 'VariableNames', {'CV_90', 'CV_95', 'CV_99'});

    fprintf('Risultati Completati.\n');
    disp(CV_Table);

end

%% ========================================================================
%% HELPER FUNCTIONS (Incluse nello stesso file per portabilità)
%% ========================================================================

function [S_m]=S_Dates_Norm(y, datevec, T, pi, m)
    % VERSIONE CORRETTA con FIX LOOP e INDEXING
    Tpi = floor(pi*T);
    vec_norm = zeros(T,1);

    % -- m=1 --
    if m==1
        Break_Dates = datevec(1,1); 
        norm1 = norm(y(Break_Dates(1,1),:))^2;
        norm2 = norm(y(Break_Dates(1,1),:)-y(T,:))^2;

        Table_norm_dates = abs([norm1, 1, Break_Dates(1,1);
                                norm2, Break_Dates(1,1)+1, T]);
        Table_norm_dates = sortrows(Table_norm_dates, 1, 'descend'); 
        A = [Table_norm_dates(:,3) - Table_norm_dates(:,2) + 1]; 
        A = cumsum(A, 1);                                         
        Table_norm_dates = [Table_norm_dates A];

        AA = A >= Tpi; 
        AA = AA & cumsum(AA,1) < 2;
        [ind_larger,~] = find(AA);
        ind_smaller = ind_larger - 1;

        if ind_smaller < 1
            Missing = Tpi;
            Start_date = Table_norm_dates(ind_larger,2) - 1;
            End_date = Table_norm_dates(ind_larger,3);
        else
            Missing = Tpi - A(ind_smaller,1);
            Start_date = Table_norm_dates(ind_larger,2) - 1;
            End_date = Table_norm_dates(ind_larger,3);
        end

        Length_reg = End_date - Start_date;

        if Length_reg == Missing 
            if ind_smaller < 1   
                pi_Dates = [Start_date+1, End_date]; 
            else 
                pi_Dates = [Table_norm_dates(1:ind_smaller,2:3); Start_date+1, End_date];
            end
        else        
            % --- FIX CRITICO 1: Loop esteso di +1 ---
            for i = 1:(Length_reg - Missing + 1)
                idx_start = Start_date + i - 1;
                % --- FIX CRITICO 2: Finestra esatta (-1) ---
                idx_end   = Start_date + i + Missing - 1;

                if idx_start == 0
                    vec_norm(i,1) = norm(y(idx_end,:))^2;
                else
                    vec_norm(i,1) = norm(y(idx_end,:) - y(idx_start,:))^2;
                end
            end

            [~, ind] = max(vec_norm);
            % Selezione finale con finestra corretta
            pi_Dates = [Table_norm_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1]; 
            pi_Dates = sortrows(pi_Dates, 1, 'ascend');
        end
        S_m = pi_Dates;

    % -- m=2, 3, 4, 5 --
    % La logica è identica a m=1, replicata per brevità con i break multipli.
    % Per compattezza, qui includo la logica generale che si applica se m>1
    % usando la stessa struttura corretta.
    elseif m >= 2 && m <= 5
        Break_Dates = datevec(1:m, m);

        % Costruzione tabella norme per i (m+1) regimi
        % Regime 1
        norms = [];
        regime_starts = [1; Break_Dates + 1];
        regime_ends   = [Break_Dates; T];

        for r = 1:(m+1)
            start_r = regime_starts(r);
            end_r   = regime_ends(r);
            if start_r == 1
                val = norm(y(end_r,:))^2;
            else
                val = norm(y(end_r,:) - y(start_r-1,:))^2;
            end
            norms = [norms; val, start_r, end_r];
        end

        Table_norm_dates = abs(norms);
        Table_norm_dates = sortrows(Table_norm_dates, 1, 'descend');

        A = [Table_norm_dates(:,3) - Table_norm_dates(:,2) + 1];
        A = cumsum(A, 1);
        Table_norm_dates = [Table_norm_dates A];

        AA = A >= Tpi;
        AA = AA & cumsum(AA,1) < 2;
        [ind_larger,~] = find(AA);
        ind_smaller = ind_larger - 1;

        if ind_smaller < 1
            Missing = Tpi;
            Start_date = Table_norm_dates(ind_larger,2) - 1;
            End_date = Table_norm_dates(ind_larger,3);
        else
            Missing = Tpi - A(ind_smaller,1);
            Start_date = Table_norm_dates(ind_larger,2) - 1;
            End_date = Table_norm_dates(ind_larger,3);
        end

        Length_reg = End_date - Start_date;

        if Length_reg == Missing
            if ind_smaller < 1
                pi_Dates = [Start_date+1, End_date];
            else
                pi_Dates = [Table_norm_dates(1:ind_smaller,2:3); Start_date+1, End_date];
            end
        else
            % --- FIX CRITICO (Applicato anche qui) ---
            vec_norm_sub = zeros(Length_reg - Missing + 1, 1);
            for i = 1:(Length_reg - Missing + 1)
                idx_start = Start_date + i - 1;
                idx_end   = Start_date + i + Missing - 1;

                if idx_start == 0
                    vec_norm_sub(i,1) = norm(y(idx_end,:))^2;
                else
                    vec_norm_sub(i,1) = norm(y(idx_end,:) - y(idx_start,:))^2;
                end
            end
            [~, ind] = max(vec_norm_sub);
            pi_Dates = [Table_norm_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
            pi_Dates = sortrows(pi_Dates, 1, 'ascend');
        end
        S_m = pi_Dates;
    end
end

function [S_m_New] = merge_Dates(S_m)
    S_m_New = [];
    i = 1;
    while i <= size(S_m,1)
        if S_m(i,1) == 1
            if S_m(i,2) ~= S_m(i+1,1) && S_m(i,2) ~= S_m(i+1,1) - 1
                S_m_New = [S_m_New; S_m(i,1), S_m(i,2)];
                i = i + 1;
            else
                j = i;
                kk = 0;
                while (S_m(j,2) == S_m(j+1,1) || S_m(j,2) == S_m(j+1,1) - 1) && kk == 0 
                    j = j + 1;
                    if j == size(S_m,1)
                        kk = j;
                        j = 1;  
                    end
                end
                if kk ~= 0
                    j = kk;
                end
                S_m_New = [S_m_New; S_m(i,1), S_m(j,2)];
                i = j + 1;
            end
        else
            if i == size(S_m,1)
                S_m_New = [S_m_New; S_m(i,1), S_m(i,2)];
                i = i + 1;
            else
                if S_m(i,2) ~= S_m(i+1,1) && S_m(i,2) ~= S_m(i+1,1) - 1
                    S_m_New = [S_m_New; S_m(i,1), S_m(i,2)];
                    i = i + 1;
                else
                    j = i;
                    kk = 0;
                    while j < size(S_m,1) && (S_m(j,2) == S_m(j+1,1) || S_m(j,2) == S_m(j+1,1) - 1) && kk == 0
                        j = j + 1;
                    end
                    S_m_New = [S_m_New; S_m(i,1), S_m(j,2)];
                    i = j + 1;
                end
            end
        end
    end
end

