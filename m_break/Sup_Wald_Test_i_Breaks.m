function SupFtest = Sup_Wald_Test_i_Breaks(y, z, i, q, bigt, datevec, prewhit, robust, x, p, hetdat, hetvar)

% NAME            : Sup_Wald_Test_i_Breaks.m
% Description     : Compute Sup_Wald_Test for fixed alternative with i breaks 
% Usage           : Sup_Wald_Test_i_Breaks(y, z, i, q, bigt, datevec, prewhit, robust, x, p, hetdat, hetvar)
% Input           : y - dependent variable                                       
% Input           : z - matrix of regressors of dimension  Txq
% Input           : i - # of breaks under alternative
% Input           : q - # of regressors
% Input           : bigt - sample size
% Input           : datevec - Break dates obtained from Global method
% Input           : prewhit - pre-whitening                                    
% Input           : robust -  Robust method for heterogeneity and autocorrelation 
%                   in the residuals, 0 otherwise. The method used is Andrews(1991) 
%                   automatic bandwidth with AR(1) approximation and the quadratic 
%                   kernel. 
% Input           : x - regressors for partial change model
% Input           : p - # of regressors in x
% Input           : hetdat - hetdat=1; % option for the construction of the F tests. Set to 1 if want to
%                    allow different moment matrices of the regressors across segments. 
%                    If hetdat=0, the same moment matrices are assumed for each segment 
%                    and estimated from the ful sample. 
%                    It is recommended to set hetdat=1 if p>0.  regressors for partial change model
% Input           : hetvar - hetvar=1;  option for the construction of the F tests. Set to 1 if want to allow 
%                   for the variance of the residuals to be different across segments.
%                   If hetvar=0, the variance of the residuals is assumed constant 
%                   across segments and constructed from the ful sample. This option 
%                   is not available when robust=1.  
% Output          : SupFtest
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
% procedure that computes the supF test for i breaks.
% construct the R matrix.

rsub=zeros(i,i+1);                                                              % i is the # of breaks under H1
j=1;
while j <= i
    rsub(j,j )= -1;                                                             % rsub is 1x2 = [-1,1]
    rsub(j,j+1) = 1;
    j =j+1;
end
rmat = kron(rsub,eye(q));
zbar = Partition_zbar(z, i, datevec(1:i,i));                                              % NAME            : Partition_zbar.m
                                                                               % Description     : procedure to construct the diagonal partition of z with
                                                                               %                   m breaks at dates bb. (bb is a column vector)
                                                                               % Usage           : zb =  Partition_zbar(zz, m, bb)            
                                                                               % Input           : zz - matrix to partition                                       
                                                                               % Input           : m - # of breaks
                                                                               % Input           : bb - column vector of dates of break 
                                                                               % Output          : zb - a Partitioned matrix of dimension zb = zeros(nt,(m+1)*q1), 
                                                                               %                   where [nt,q1] = size(zz);   
%==========================================================================
%            Pure SC Model
%==========================================================================                                                                               
if p==0
    delta = ols_simple(y,zbar);                                                 % Esimate by ols
else
%==========================================================================
%            Partial SC Model
%==========================================================================
    dbdel = olsqr(y,[zbar x]);
    delta = dbdel(1:(i+1)*q,1);                                                 % Recall i breask => i+1 regimes
end

vdel = Var_delta(y,z,i,q,bigt,datevec(1:i,i),prewhit,robust,x,p,0,hetdat,hetvar);

fstar = delta'*rmat'*inv(rmat*vdel*rmat')*rmat*delta;
SupFtest = (bigt-(i+1)*q-p)*fstar/(bigt*i);                                     % Wald Sup Test 

