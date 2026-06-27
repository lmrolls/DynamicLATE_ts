function [beta_hat, S_vec, dates] = est_fgls_pilate(Y, D, Z, X, pi0, m)
    % Setup del sistema multivariato
    y_bar = [Y, D]; 
    little_Z = [Z, X];
    q = size(Z, 2); 
    p = size(X, 2);
    T = size(Y, 1);
    R_selection = eye(2 * (q + p));
    
    % 1. Identificazione break multivariata (Qu & Perron)
    datevec = mbreak_multivariate(y_bar, little_Z, R_selection, m);
    [maty, matz] = transf(y_bar, little_Z);
    T_multi = size(maty, 1);
    
    % 2. Selezione Greedy multivariata
    dates_multi = S_Dates_multi_greedy_double_eq(maty, matz, datevec, T_multi, pi0, m, q, p);
    
    % 3. Stima FGLS e generazione output
    [~, Xi_hat, S_index_multi] = FGLS_estimation_S(maty, matz, floor(pi0*T), q, p, dates_multi);
    
    % Estrazione beta e mappatura S
    beta_hat = extract_beta(Xi_hat, q, p);
    S_vec = S_index_multi(1:2:end);
    dates = [ (dates_multi(:,1)+1)/2, dates_multi(:,2)/2 ];
end


function beta = extract_beta(Xi, q, p)
    if p == 0, beta = Xi(1,1) / Xi(q+1,1); 
    else, beta = Xi(1,1) / Xi(p+q+1,1); end
end


function [Obj, Xi_hat_S, S_index] = FGLS_estimation_S(maty, matz_ST, T_pi, q, p, Sup_S)
    T_multi = size(maty, 1); S_index = zeros(T_multi, 1);
    for irow = 1:size(Sup_S, 1), S_index(Sup_S(irow,1):Sup_S(irow,2)) = 1; end
    matz_S = matz_ST; matz_S(S_index == 0, :) = 0;
    [Xi_hat] = OLS_Robust_multi_full(maty, matz_S, [], 0, q, p);
    epsilon_hat = maty - matz_S * Xi_hat;
    idx_odd = 1:2:T_multi; active_mask_t = S_index(idx_odd) == 1; 
    u_S = epsilon_hat(idx_odd(active_mask_t)); v_S = epsilon_hat(1+idx_odd(active_mask_t));
    E_S = [u_S, v_S]; T_pi_eff = length(u_S);
    Sigma_hat = (E_S' * E_S) / (T_pi_eff - q - p); invSigma = inv(Sigma_hat);
    w11 = invSigma(1,1); w12 = invSigma(1,2); w22 = invSigma(2,2);
    Y1 = maty(1:2:end); Y2 = maty(2:2:end); Z1 = matz_S(1:2:end, :); Z2 = matz_S(2:2:end, :);
    Y1(~active_mask_t) = 0; Y2(~active_mask_t) = 0;
    M_XX = w11 * (Z1'*Z1) + w12 * (Z1'*Z2 + Z2'*Z1) + w22 * (Z2'*Z2);
    M_XY = w11 * (Z1'*Y1) + w12 * (Z1'*Y2) + w12 * (Z2'*Y1) + w22 * (Z2'*Y2);
    Xi_hat_S = M_XX \ M_XY;
    res_S = maty - matz_S * Xi_hat_S; res1 = res_S(1:2:end); res2 = res_S(2:2:end);
    res1(~active_mask_t) = 0; res2(~active_mask_t) = 0;
    Obj = w11*(res1'*res1) + 2*w12*(res1'*res2) + w22*(res2'*res2);
end

function Sup_S = S_Dates_multi_greedy_double_eq(maty, matz, datevec, T_multi, pi, m, q, p)
    % Dimensioni in scala singola (T) e numero minimo di osservazioni (Tpi)
    T = floor(T_multi / 2); 
    Tpi = floor(pi * T);
    
    % --- FIX CRITICO: Estrazione corretta della colonna m da datevec ---
    if isempty(datevec) || m == 0
        datevec = []; 
    else
        % Prende SOLO le date del modello con m break (evita la frammentazione)
        datevec = datevec(1:m, m); 
        datevec = sort(unique(datevec(datevec>0 & datevec<T))); 
    end
    
    Full_B = [0; datevec; T]; 
    Num_R = length(Full_B) - 1;
    Table_R = zeros(Num_R, 4);
    
    % Calcolo Statistica F_sys per ogni regime
    for i = 1:Num_R
        S = Full_B(i)+1; 
        E = Full_B(i+1); 
        L = E-S+1; 
        
        % Mappatura esatta sulla scala 2T per i dati interleaved
        idx_2T = (S-1)*2+1 : E*2;
        [~, ~, F_sys] = OLS_Robust_multi_full(maty(idx_2T), matz(idx_2T,:), [], 0, q, p);
        Table_R(i,:) = [S, E, L, F_sys];
    end
    
    % Ordinamento Greedy: dal regime col segnale più forte al più debole
    Table_R = sortrows(Table_R, 4, 'descend');
    A = cumsum(Table_R(:,3)); 
    idx_L = find(A >= Tpi, 1, 'first');
    
    if isempty(idx_L)
        idx_L = size(Table_R, 1); 
    end
    
    idx_S = idx_L - 1; 
    Missing = Tpi; 
    if idx_S >= 1
        Missing = Tpi - A(idx_S); 
    end
    
    S_t = Table_R(idx_L,1); 
    E_t = Table_R(idx_L,2); 
    L_t = Table_R(idx_L,3);
    
    % Selezione Finale / Sliding Window
    if L_t == Missing
        % Incastro perfetto
        if idx_S < 1
            pi_D = Table_R(idx_L, 1:2); 
        else
            pi_D = Table_R(1:idx_L, 1:2); 
        end
    else
        % Ricerca della finestra ottima (Sliding Window senza errori di indice)
        max_f = -1; 
        b_off = 0;
        
        for k = 0:(L_t-Missing)
            % Indici precisi per estrarre coppie intatte Y/D nella scala 2T
            sub_idx = (S_t+k-1)*2+1 : (S_t+k+Missing-1)*2;
            [~,~,f] = OLS_Robust_multi_full(maty(sub_idx), matz(sub_idx,:), [], 0, q, p);
            
            if f > max_f
                max_f = f; 
                b_off = k; 
            end
        end
        
        Best_W = [S_t+b_off, S_t+b_off+Missing-1];
        
        if idx_S < 1
            pi_D = Best_W; 
        else
            pi_D = [Table_R(1:idx_S, 1:2); Best_W]; 
        end
    end
    
    % Riordino cronologico e riconversione in scala 2T per l'output SUR/FGLS
    pi_D = sortrows(pi_D, 1); 
    Sup_S = [(pi_D(:,1)-1)*2+1, pi_D(:,2)*2];
end


function [estimator, t_stat, Fstat] = OLS_Robust_multi_full(Y_s, Z_s, X_s, HAC, q, p)
    [n_2T, nvar] = size(Z_s); T = n_2T/2;
    XX = Z_s' * Z_s;
    if rank(XX) < nvar, estimator = NaN(nvar,1); t_stat = 0; Fstat = 0; return; end
    estimator = XX \ (Z_s' * Y_s);
    res = Y_s - Z_s * estimator;
    meat = zeros(nvar, nvar);
    for t = 1:T
        idx = (t-1)*2 + (1:2); Z_t = Z_s(idx,:); r_t = res(idx);
        meat = meat + Z_t' * (r_t * r_t') * Z_t;
    end
    vcov = (XX\eye(nvar)) * meat * (XX\eye(nvar)) * (n_2T/(n_2T-nvar));
    t_stat = abs(estimator ./ sqrt(diag(vcov)));
    idx_theta = (p+q+1) : (p+q+q); Fstat = mean(t_stat(idx_theta).^2); 
end