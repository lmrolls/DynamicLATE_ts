function [S_m] = S_Dates_m_select(y, z, x, T, Tpi, Table, A)
% S_DATES_M_SELECT Sub-sample of Tpi dates from the sorted regime table (step 2 of S_Dates_m).
%
%   S_m = S_Dates_m_select(y, z, x, T, Tpi, Table, A)
%
% Takes the regimes in decreasing order of |stat| (Table from S_Dates_m_regimes) until the
% cumulative length reaches Tpi; the last regime is entered only partially, through the window
% of 'Missing' dates with the largest |t| (q = 1) or F (q > 1) statistic among the sliding
% windows inside that regime. The window loop is the one of the original S_Dates_m (branch
% tryMM): windows y(Start_date+i : Start_date+i+Missing), i = 1, ..., Length_reg - Missing.
%
% OUTPUT: S_m (n x 2) = [start, end] of the selected regimes, sorted by start date.

    q = size(z, 2);
    vec_tstat = zeros(T, 1);
    vec_Fstat = zeros(T, 1);

    AA = A >= Tpi;
    AA = AA & cumsum(AA, 1) < 2;
    [ind_larger, ~] = find(AA);
    ind_smaller = ind_larger - 1;
    if ind_smaller < 1
        Missing = Tpi;
        Start_date = Table(ind_larger, 2) - 1;
        End_date = Table(ind_larger, 3);
    else
        Missing = Tpi - A(ind_smaller, 1);
        Start_date = Table(ind_larger, 2) - 1;
        End_date = Table(ind_larger, 3);
    end
    Length_reg = End_date - Start_date;
    if Length_reg == Missing
        if ind_smaller < 1
            pi_Dates = [Start_date + 1, End_date];
        else
            pi_Dates = [Table(1:ind_smaller, 2:3); Start_date + 1, End_date];
            pi_Dates = sortrows(pi_Dates, 1, 'ascend');
        end
    else
        for i = 1:Length_reg - Missing
            idx = Start_date + i:Start_date + i + Missing;
            if isempty(x) == 1
                if q == 1
                    [~, vec_tstat(i, 1)] = OLS_Robust(y(idx), z(idx, :), [], 0);
                else
                    [~, ~, vec_Fstat(i, 1)] = OLS_Robust(y(idx), z(idx, :), [], 0);
                end
            else
                if q == 1
                    [~, vec_tstat(i, 1)] = OLS_Robust(y(idx), z(idx, :), x(idx, :), 0);
                else
                    [~, ~, vec_Fstat(i, 1)] = OLS_Robust(y(idx), z(idx, :), x(idx, :), 0);
                end
            end
        end
        if q == 1
            [~, ind] = max(abs(vec_tstat));
        else
            [~, ind] = max(abs(vec_Fstat));
        end
        pi_Dates = [Table(1:ind_smaller, 2:3); Start_date + ind, Start_date + ind + Missing - 1];
        pi_Dates = sortrows(pi_Dates, 1, 'ascend');
    end
    S_m = pi_Dates;
end
