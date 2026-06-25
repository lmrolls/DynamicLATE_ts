function [maxf, newd]= SeqTest(bigvec, dt, nseg, y, z, h, q, prewhit, robust, x, p, hetdat, hetvar)

% NAME            : SeqTest.m
% Description     : procedure that computes the supF(l+1|l) test l=nseg-1. The l breaks under
%                   the null are taken from the global minimization (in dt).
% Usage           : Sup_Wald_Test_i_Breaks(y, z, i, q, bigt, datevec, prewhit, robust, x, p, hetdat, hetvar)
% Input           : bigvec - big vector of recursive residuals 
% Input           : dt - break dates from glaobal method
% Input           : y - dependent variable                                       
% Input           : z - matrix of regressors of dimension  Txq
% Input           : h - minimum regime length
% Input           : q - # of regressors
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
% Output          : maxf  - SupFtest for testing l+1|l breaks
% Output          : newd - new break date
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
% procedure that computes the supF(l+1|l) test l=nseg-1. The l breaks under
% the null are taken from the global minimization (in dt).

ssr = zeros(nseg,1);                                                              % nseg is # of regimes 
ftestv = zeros(nseg,1);
bigt = size(z,1);                                                               % sample size
dv=zeros(nseg+1,1);
dv(2:nseg,1) = dt;
dv(nseg+1,1) = bigt;
ds = zeros(nseg,1);

in=0;

for is = 1:nseg
    length = dv(is+1,1)-dv(is,1);                                               % size of isth regime
    
    if length >= 2*h
        if p==0
            [ssr(is,1), ds(is,1)] = Partition_1_Break(dv(is,1)+1,...
                dv(is,1)+h, dv(is+1,1)-h, dv(is+1,1),bigvec, bigt);                % NAME            : Partition_1_Break.m
                                                                                % Description     : Obtain an optimal one break partition (using Recursive Residual) for a segment that
                                                                                %                   starts at date "start" and ends at date "last". It returns the optimal break
                                                                                %                   and the associated SSR.     
                                                                                % Usage           : Partition_1 (start,b1,b2,last,bigvec,bigt)
                                                                                % Input           : start - beginning of the segment considered
                                                                                % Input           : last - ending date of the last segment considered.
                                                                                % Input           : b1 - first possible break date (usually set to h)                                    
                                                                                % Input           : b2: last possible break date
                                                                                % Input           : bigvec - vector storing recursive residuals (output
                                                                                %                   from Recursive_SSR.m
                                                                                % Input           : bigt - Sample size
                                                                                % Output          : ssrmin - Minimum sum of squared residuals (SSR) 
                                                                                % Output          : dx - Optimal break
                                                                                
                                                                                % dv(is,1)+1 = start date
                                                                                % dv(is+1,1) = last
                                                                                % b1 - first possible break date (usually set to h) is at dv(is,1)+1
                                                                                % b2: last possible break date: dv(is+1,1)-h
                                                                                
            ftestv(is,1)= Sup_Wald_Test_i_Breaks(y(dv(is,1)+1:dv(is+1,1),1), z(dv(is,1)+1:dv(is+1,1),:),...
                1, q, length, ds(is,1)-dv(is,1), prewhit,robust,0,p, hetdat, hetvar);
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
            
        else
            [ssr(is,1), ds(is,1)] = onebp(y,z,x,h,dv(is,1)+1,dv(is+1,1));
            ftestv(is,1) = Sup_Wald_Test_i_Breaks(y(dv(is,1)+1:dv(is+1,1),1), z(dv(is,1)+1:dv(is+1,1),:),...
                1,q,length,ds(is,1)-dv(is,1),prewhit,robust,x(dv(is,1)+1:dv(is+1,1),:),...
                p,hetdat,hetvar);
        end
        
    else
        in=in+1;
        ftestv(is,1)=0.0;
    end
end

if in == nseg
    disp(['Given the location of the breaks from the global optimization with ' ...
        num2str(1) ' breaks there was no more place to insert' ...
        ' an additional breaks that satisfy the minimal length requirement.']);
end

[maxf,news] = max(ftestv(1:nseg,1));                                            % So I guess the l+1 is the one whihc maximizes ftestv
newd = ds(news,1);



        
            