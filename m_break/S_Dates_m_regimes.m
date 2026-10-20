function [Table, A] = S_Dates_m_regimes(y, z, x, datevec, T, m)
% S_DATES_M_REGIMES Regime-level statistics of the sub-sample selection rule (step 1 of S_Dates_m).
%
%   [Table, A] = S_Dates_m_regimes(y, z, x, datevec, T, m)
%
% Splits {1,...,T} into the m+1 regimes delimited by the Bai-Perron break dates datevec(1:m, m),
% computes on each regime the robust t-statistic of z (q = 1) or the F-statistic (q > 1) with
% OLS_Robust, and returns the table [|stat|, start, end, cumulative length] sorted by |stat| in
% descending order, exactly as the m = 1, ..., 5 blocks of the original S_Dates_m did.
%
% This table depends on (y, z, x, datevec, m) only, not on pi: S_Dates_m_grid computes it once
% per m and reuses it for every pi of the search grid (same operations, same order, identical
% results). The selection for a given Tpi is done by S_Dates_m_select.
%
% INPUTS : y (T x 1), z (T x q), x (T x p or []), datevec (triangular matrix of break dates as
%          returned by Procedure_piL_Full_Sample), T, m (number of breaks, 1 <= m <= size(datevec,2))
% OUTPUTS: Table (m+1 x 4) = [|stat|, start, end, cumulative length], sorted; A = Table(:,4)

    bd = datevec(1:m, m);
    starts = [1; bd(:) + 1];
    ends   = [bd(:); T];
    q = size(z, 2);
    stat = zeros(m + 1, 1);
    for i = 1:m + 1
        idx = starts(i):ends(i);
        if isempty(x) == 1
            [~, t_stat_i, Fstat_i] = OLS_Robust(y(idx), z(idx, :), [], 0);
        else
            [~, t_stat_i, Fstat_i] = OLS_Robust(y(idx), z(idx, :), x(idx, :), 0);
        end
        if q == 1
            stat(i) = t_stat_i;
        else
            stat(i) = Fstat_i;
        end
    end
    Table = abs([stat, starts, ends]);
    Table = sortrows(Table, 1, 'descend');
    A = [Table(:, 3) - Table(:, 2) + 1];
    A = cumsum(A, 1);
    Table = [Table A];
end
