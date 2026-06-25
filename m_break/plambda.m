function lambda= Partitioned_lambda(b, m, bigt)

% NAME            : Partitioned_lambda.m
% Description     : Construct a diagonal matrix of dimension m+1 with ith
%                   entry (T_i-T_i-1)/T.
% Usage           : Var_delat(y, z, i, q, bigt, b, prewhit, robust, x, p, withb, hetdat, hetvar)
% Input           : m - # of breaks under alternative 
% Input           : bigt - sample size
% Input           : b - Break dates obtained from Global method
% Output          : lambda
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
% procedure that compute the covariance matrix of the estimates delta.
% procedure that construct a diagonal matrix of dimension m+1 with ith
% entry (T_i-T_i-1)/T.

lambda = zeros(m+1,m+1);                                                        % m+1 regimes
lambda(1,1) = b(1,1)/bigt;
k=2;
while k <= m
    lambda(k,k)=(b(k,1)-b(k-1,1))/bigt;
    k=k+1;
end
lambda(m+1,m+1) = (bigt-b(m,1))/bigt;
