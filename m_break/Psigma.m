function sigmat = Psigma(res,b,q,i,nt)

% NAME            : Psigma.m
% Description     : Compute estimate of variance-cov of residuals while it
%                   is allowed to change across regimes. It returns a
%                   (i+1) diagonal matrix with ith entry the estimate of 
                    % the variance of residuals in the ith regime. 
% Usage           : sigmat = psigmq(res,b,q,i,nt)
% Input           : res - residuals                                       
% Input           : b - Break dates obtained from Global method
% Input           : q - # of regressors
% Input           : i - # of breaks under alternative
% Input           : nt - sample size
% Output          : sigmat
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

sigmat = zeros(i+1,i+1);
sigmat(1,1)=res(1:b(1,1),1)'*res(1:b(1,1),1)/b(1,1);
kk = 2;
while kk <= i
    sigmat(kk,kk) = res(b(kk-1,1)+1:b(kk,1),1)'*res(b(kk-1,1)+1:b(kk,1),1)...
        /(b(kk,1)-b(kk-1,1));                                                   % Because if b1=5 and b2=10 then 10-5
    kk=kk+1;
end

sigmat(i+1,i+1)=res(b(i,1)+1:nt,1)'*res(b(i,1)+1:nt,1)/(nt-b(i,1));             % Last regime
