function [ssrmin, dx]= Partition_1_Break(start, b1, b2, last, bigvec, bigt)

% NAME            : Partition_1_Break.m
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
% procedure to obtain an optimal one break partitions for a segment that
% starts at date start and ends at date last. It returns the optimal break
% and the associated SSR.

% start: beginning of the segment considered
% b1: first possible break date 
% b2: last possible break date
% last: end of segment considered

dvec = zeros(bigt,1);
ini = (start-1)*bigt - (start-2)*(start-1)/2 + 1;
j=b1;
while j <= b2
    jj= j - start;
    k=j*bigt - (j-1)*j/2 + last-j;
    dvec(j,1)=bigvec(ini+jj,1)+ bigvec(k,1);
    j=j+1;
end

ssrmin = min(dvec(b1:b2,1))';                                                   % ssrmin is minimum srr for one break
[dummy, minindcdvec]=min(dvec(b1:b2,1));                                        % [C,I] = min(...) finds the indices of the minimum values of A, and returns them in output vector I.
                                                                                % If there are several identical minimum values, the index of the first one found is returned.
                                                                                % minindcdvec is index of minimum
dx = (b1-1) + minindcdvec';

    
    