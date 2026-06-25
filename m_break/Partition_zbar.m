function zb = Partition_zbar(zz, m, bb)

% NAME            : Partition_zbar.m
% Description     : procedure to construct the diagonal partition of z with
%                   m breaks at dates bb. (bb is a column vector)
% Usage           : zb =  Partition_zbar(zz, m, bb)            
% Input           : zz - matrix to partition                                       
% Input           : m - # of breaks
% Input           : bb - column vector of dates of break 
% Output          : zb - a Partitioned matrix of dimension zb = zeros(nt,(m+1)*q1), 
%                   where [nt,q1] = size(zz);          
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

% procedure to construct the diagonal partition of z with m break at date
% b.

[nt,q1] = size(zz);                                                             % size of regressors

zb = zeros(nt,(m+1)*q1);                                                        % zbar is zb. zb is ntx(m+1)*q1 (because if there are m breaks than there are m+1 regimes                                                  
zb(1:bb(1,1),1:q1)=zz(1:bb(1,1),:);                                             % Fix first paritition
i=2;                                                                            % fix 2nd breaks and do loop with while
while i <= m
    zb(bb(i-1,1)+1:bb(i,1),(i-1)*q1+1:i*q1) = zz(bb(i-1,1)+1:bb(i,1),:);
    i=i+1;
end
zb(bb(m,1)+1:nt, m*q1+1:(m+1)*q1) = zz(bb(m,1) + 1:nt,:);                       % Store last regime
