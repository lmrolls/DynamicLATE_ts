function [ssrmin,dx] = Partition_1_Seq(b1, b2, last, vssr, vssrev)

% NAME            : Partition_1_Seq.m
% Description     : procedure to obtain an optimal one break partitions for a segment that
%                   starts at date start and ends at date last. It returns the optimal break
%                   and the associated SSR.
%                   Procedure used with the sequential method whenthe T*(T+1)/2 vector of SSR
%                   is not computed.      
% Usage           : Partition_1_Seq(b1,b2,last,vssr,vssrev)
% Input           : start - starting entry of the sample used.
% Input           : last - ending date of the last segment considered.
% Input           : b1 - first possible break date (usually set to h)                                    
% Input           : b2: last possible break date
% Input           : vssr - vector storing recursive residuals (output
%                   from Sequential_Method.m). These are the pre-break SSR.
% Input           : vssrev - vector storing recursive residuals (output
%                   from Sequential_Method.m). These are the pre-break SSR.
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
% this procedure computes recursive residuals from a data set that starts
% at date "start" and ends at date "last". It returns a vector of sum of
% squared residuals (SSR) of length last-start+1 (stored for convenience in
% a vector of length T).

% start: starting entry of the sample used.
% last: ending date of the last segment considered.
% y: dependent variable
% z: matrix of regressors of dimension Txq
% h: minimal length of a segment


% procedure to obtain an optimal one break partitions for a segment that
% starts at date start and ends at date last. It returns the optimal break
% and the associated SSR.
% Procedure used with the sequential method whenthe T*(T+1)/2 vector of SSR
% is not computed.

% start: beginning of the segment considered
% b1: first possible break date
% b2: last possible break date
% last end of segment considered

dvec = zeros(last,1);

j=b1;

while j <= b2
    dvec(j,1)= vssr(j,1)+ vssrev(last-j,1);                                     % SSR from pre-break and post break SSR
    j=j+1;
end

[ssrmin ssrminind] = min(dvec(b1:b2,1));                                        % Pick the minimum
dx=(b1-1) + ssrminind;                                                          % Because if break at b1 then the index is 1 and so b1-1+1=b1
