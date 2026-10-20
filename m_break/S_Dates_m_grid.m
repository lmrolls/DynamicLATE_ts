function [S_cells, Tpi_list] = S_Dates_m_grid(y, z, x, datevec, T, Tpi_list, m_list)
% S_DATES_M_GRID Sub-sample selection for a whole grid of (Tpi, m) without repeated work.
%
%   [S_cells, Tpi_list] = S_Dates_m_grid(y, z, x, datevec, T, Tpi_list, m_list)
%
% S_cells{k, j} is the matrix of selected regimes for Tpi = Tpi_list(k) dates and m = m_list(j)
% breaks, identical to S_Dates_m(y, z, x, datevec, T, pi, m) with floor(pi*T) = Tpi_list(k).
%
% Computational savings with respect to calling S_Dates_m for every cell of the search grid of
% piLATE_Ftest / est_2sls_pilate (pi_low x pi x m):
%   - the regime-level statistics (m+1 regressions) are computed once per m instead of once per
%     (pi, m) cell: S_Dates_m_regimes;
%   - each distinct Tpi is processed once per m, whatever the number of pi values of the grid
%     that map to it: S_Dates_m_select.
% Every statistic is computed with the same function calls, on the same data, in the same order
% as in S_Dates_m, so the results coincide bit for bit.
%
% INPUTS : y (T x 1), z (T x q), x (T x p or []), datevec (from Procedure_piL_Full_Sample), T,
%          Tpi_list (vector of integers, e.g. floor(pi_grid * T)), m_list (vector of integers >= 1)
% OUTPUT : S_cells (numel(Tpi_list) x numel(m_list) cell array)

    nP = numel(Tpi_list);
    nM = numel(m_list);
    S_cells = cell(nP, nM);
    for j = 1:nM
        m = m_list(j);
        [Table, A] = S_Dates_m_regimes(y, z, x, datevec, T, m);
        done_Tpi = [];
        done_S = {};
        for k = 1:nP
            Tpi = Tpi_list(k);
            kk = find(done_Tpi == Tpi, 1);
            if isempty(kk)
                S_cells{k, j} = S_Dates_m_select(y, z, x, T, Tpi, Table, A);
                done_Tpi(end + 1) = Tpi; %#ok<AGROW>
                done_S{end + 1} = S_cells{k, j}; %#ok<AGROW>
            else
                S_cells{k, j} = done_S{kk};
            end
        end
    end
end
