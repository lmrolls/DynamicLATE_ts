function Results = piLATE_Ftest(D, Z, X, options)
% piLATE_Ftest Executes the pi-LATE (F*) test for structural identification.
%
% SYNTAX:
%   Results = piLATE_Ftest(D, Z, [], options)  % Without exogenous
%   covariates
%   Results = piLATE_Ftest(D, Z, X, options)   % With exogenous covariates
%
% INPUTS:
%   D       : Vector of the endogenous variable (T x 1)
%   Z       : Matrix of instruments (T x q)
%   X       : Matrix of exogenous variables (T x p). Can be empty [].
%   options : Optional struct with algorithm parameters:
%             .piLow              (Default: 0.6)
%             .Hac_option         (Default: false)
%             .m_breaks           (Default: 5)
%
% OUTPUT:
%   Results : Struct containing:
%             .F_star       : pi-LATE F* statistic
%             .best_pi      : Optimal sample fraction
%             .best_m       : Optimal number of breaks
%             .F_FullSample : F-statistic on the full sample
%             .S_index_opt  : Logical vector (T x 1) indicating the sub-sample
%
% =========================================================================
    warning('off','all');
    % --- OPTIONS AND DEFAULT HANDLING ---
    if nargin < 3, X = []; end
    if nargin < 4 || isempty(options), options = struct(); end
    if ~isfield(options, 'piLow'),              options.piLow = 0.6; end
    if ~isfield(options, 'Hac_option'),         options.Hac_option = false; end
    if ~isfield(options, 'm_breaks'),           options.m_breaks = 5; end
    
    % Fixed Bai-Perron parameters (paper standard)
    eps1 = 0.05; fixb = 0; betaini = 0; maxi = 20; printd = 0; eps_val = 0.0001; 
    robust = 1; prewhit = 0; hetdat = 1; hetvar = 1; hetomega = 1; hetq = 1; doglobal = 1;
    T = size(D, 1);
    q = size(Z, 2);
    if isempty(X), p = 0; else, p = size(X, 2); end
    h = round(eps1 * T);
    
    % --- 1. BREAKS IDENTIFICATION (Bai-Perron / External dependency) ---
    [~, datevec, ~] = Procedure_piL_Full_Sample(T, D, Z, q, options.m_breaks, h, eps1, ...
        robust, prewhit, hetomega, hetq, doglobal, fixb, X, p, eps_val, maxi, betaini, printd, hetdat, hetvar);
        
    % --- 2. SEARCH FOR THE F* STATISTIC ---
    vector_piL = options.piLow : 0.1 : 1;
    step_pi = 0.05;
    
    Sup_F_piL  = zeros(length(vector_piL), 1);
    Sup_m_piL  = zeros(length(vector_piL), 1);
    Sup_pi_piL = zeros(length(vector_piL), 1);
    for piL_idx = 1:length(vector_piL)
        num_pis = vector_piL(piL_idx) : step_pi : 1;
        F_ST = zeros(length(num_pis), options.m_breaks);
        
        for m = 1:options.m_breaks
            for pi_idx = 1:length(num_pis)
                % External dependency
                S_m = S_Dates_m(D, Z, X, datevec, T, num_pis(pi_idx), m);
                S_index = zeros(T, 1);
                for irow = 1:size(S_m, 1), S_index(S_m(irow,1):S_m(irow,2)) = 1; end
                
                % Compression
                [D_tilde, Z_tilde, T_eff] = apply_Ct_operator(D, Z, X, S_index);
                
               
                F_ST(pi_idx, m) = F_ST_Stat_Legacy(D_tilde, Z_tilde, T_eff, q, p, options.Hac_option);
                
            end
        end
        
        [Val_m, ind_m] = max(F_ST, [], 1); 
        [~, ind_m_pi] = max(Val_m); 
        Sup_m_piL(piL_idx) = ind_m_pi;
        Sup_pi_piL(piL_idx) = num_pis(ind_m(ind_m_pi));
        Sup_F_piL(piL_idx) = max(F_ST, [], 'all'); 
    end
    
    [F_star, index_Sup_piL] = max(Sup_F_piL);
    best_pi = Sup_pi_piL(index_Sup_piL);
    best_m  = Sup_m_piL(index_Sup_piL);
    
    % --- 3. FULL-SAMPLE F-STAT COMPUTATION ---
    if isempty(X)
        F_FullSample = F_ST_Stat_Legacy(D, Z, T, q, p, options.Hac_option);
    else
        % Brutal inversion for exact match with legacy floating-point precision
        Mx = eye(T) - X * (X' * X)^-1 * X';
        F_FullSample = F_ST_Stat_Legacy(Mx*D, Mx*Z, T, q, p, options.Hac_option);
    end
    
    % --- 4. OUTPUT PREPARATION ---
    S_m_opt = S_Dates_m(D, Z, X, datevec, T, best_pi, best_m);
    S_index_opt = zeros(T, 1);
    for irow = 1:size(S_m_opt, 1), S_index_opt(S_m_opt(irow,1):S_m_opt(irow,2)) = 1; end
    
    Results.F_star = F_star;
    Results.best_pi = best_pi;
    Results.best_m = best_m;
    Results.F_FullSample = F_FullSample;
    Results.S_index_opt = S_index_opt;
end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================
function [D_tilde, Z_tilde, T_eff] = apply_Ct_operator(D, Z, X, S_input)
    T = size(D, 1);
    S_index = zeros(T, 1);
    
    if isempty(S_input)
    elseif size(S_input, 2) == 2 && all(S_input(:) >= 1) && max(S_input(:)) <= T+1
        for irow = 1:size(S_input, 1), S_index(round(S_input(irow, 1)) : round(S_input(irow, 2))) = 1; end
    elseif length(S_input) == T && all(ismember(S_input(:), [0, 1]))
        S_index = S_input(:);
    elseif islogical(S_input) && length(S_input) == T
        S_index(S_input) = 1;
    else
        idx = round(S_input(:)); idx = idx(idx >= 1 & idx <= T); S_index(idx) = 1;
    end
    
    idx_active = (S_index == 1);
    T_eff = sum(idx_active);
    if T_eff == 0, D_tilde = zeros(0, size(D, 2)); Z_tilde = zeros(0, size(Z, 2)); return; end
    
    D_sub = D(idx_active, :);
    Z_sub = Z(idx_active, :);
    
    if ~isempty(X)
        X_sub = X(idx_active, :);
        % Explicit inversion ^-1 for perfect Floating Point reproduction
        Mx_sub = eye(T_eff) - X_sub * (X_sub' * X_sub)^-1 * X_sub';
        D_tilde = Mx_sub * D_sub;
        Z_tilde = Mx_sub * Z_sub;
    else
        D_tilde = D_sub;
        Z_tilde = Z_sub;
    end
end

function [F_ST_Stat] = F_ST_Stat_Legacy(D, Z, T_val, q, p, hac)
    % Classical F-stat with sample divisor
    if nargin < 6; hac = 1; end
    if T_val <= (p + q), F_ST_Stat = 0; return; end
    J_hat_ST = OLS_Robust_VarXe(D, Z, 0, hac, 0); 
    F_ST_Stat = (D'*Z) * (J_hat_ST)^-1 * (Z'*D) / (q*(T_val - p - q));
end
