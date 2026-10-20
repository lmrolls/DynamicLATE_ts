% NAME            : Run_pi_LATE.m
% Description     :        
% Usage           : 
% Input           : 
% Input           : 
% Input           : 
%                                       
% Input           : sqrt(evar_m(JI)) - square root of right window spot variance estimate
% Input           : Test_Val - is value of test stastic
% Output          : beta_eff  = Jump Beta Efficient Estimator 
%==========================================================================
% AUTHOR          : ALESSANDRO CASINI
% Email           : acasini@bu.edu  
% Date            : First version: 8/28/2015
% Last update     : 10/2/2015

%==========================================================================
% DISCLAIMER       
%==========================================================================  
% Please let me know if you find any bugs or encounter any problems while 
% using this code. All feedback is appreciated.
% @copyright, 2015

%==========================================================================
% REFERENCES
%==========================================================================
%  Bai, Jushan and Pierre Perron (1998): "Estimating and Testing Linear
%  Models with Multiple Structural Changes," Econometrica, vol 66, 
% and 
%  Bai, Jushan and Pierre Perron (2003): "Computation and Analysis of
%  Multiple Structural Change Models," Journal of Applied Econometrics, 18, 1-22.
% 
%
%==========================================================================
% DESCRIPTIONS: Computation Method
%==========================================================================
%
%
%==========================================================================
% COMMENTS
%==========================================================================

%% ========================================================================
%             Preliminaries
%  ========================================================================
% close all; clear all; clc;

function [S_m]=S_Dates_m(y, z, x, datevec, T, pi,m)
% S_DATES_M Sub-sample of floor(pi*T) dates selected from the m+1 Bai-Perron regimes.
%
% [branch perf] The five copy-pasted blocks for m = 1, ..., 5 of the original file are replaced
% by one general implementation in two steps, valid for any m <= size(datevec, 2):
%   S_Dates_m_regimes : regime-level |t| (q = 1) or F (q > 1) statistics, sorted table;
%   S_Dates_m_select  : greedy selection of regimes up to Tpi dates, sliding window inside the
%                       last (partially included) regime.
% The operations, their order and the window convention are those of the original blocks, so the
% output is identical; S_Dates_m_grid uses the same two steps for a whole grid of (Tpi, m)
% without recomputing the regime statistics for every pi.

Tpi=floor(pi*T);
[Table, A] = S_Dates_m_regimes(y, z, x, datevec, T, m);
S_m = S_Dates_m_select(y, z, x, T, Tpi, Table, A);
end