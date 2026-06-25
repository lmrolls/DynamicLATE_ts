% NAME            : Run_pi_LATE.m
% Description     :        
% Usage           : 
% Input           : 
% Input           : 
% Input           : 
%                                       
% Input           : 
% Input           : 
% Output          : 
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
% DESCRIPTIONS: 
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

function [F_ST_Stat]=F_ST_Stat(D, Z, Tpi, q, p, hac)
if nargin < 6;    hac=1; end
%J_hat_ST=OLS_Robust_VarXe(D,Z,0, 0, 0);                                     % Heteroskedasticity-only
J_hat_ST=OLS_Robust_VarXe(D,Z,0, hac, 0);                                    % Hac
F_ST_Stat=D'*Z*(J_hat_ST)^-1 * Z'*D /(q*(Tpi - p- q));
end