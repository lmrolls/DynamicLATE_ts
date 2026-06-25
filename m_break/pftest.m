function ftest = SupFTest_i_Breaks(y,z,i,q,bigt,datevec,prewhit,robust,x,p,hetdat,hetvar)

% NAME            : SupFTest_i_Breaks.m
% Description     : procedure that computes the supF test for i breaks.
%                   construct the R matrix.
% Usage           : Recursive_SSR(start,y,z,h,last)
% Input           : start - starting entry of the sample used.
% Input           : last - ending date of the last segment considered.
% Input           : y - dependent variable                                       
% Input           : z - matrix of regressors of dimension  Txq
% Input           : h - minimal length of a segment
% Output          : vec_Recursive_SSR - vector of sum of
%                   squared residuals (SSR) of length last-start+1 (stored for convenience in
%                   a vector of length T). 
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


% construct the R matrix.

rsub=zeros(i,i+1);
j=1;
while j <= i
    rsub(j,j)=-1;
    rsub(j,j+1)=1;
    j=j+1;
end
rmat=kron(rsub,eye(q));
zbar=pzbar(z,i,datevec(1:i,i));
if p==0
    delta=olsqr(y,zbar);
else
    dbdel=olsqr(y,[zbar x]);
    delta=dbdel(1:(i+1)*q,1);
end

vdel=pvdel(y,z,i,q,bigt,datevec(1:i,i),prewhit,robust,x,p,0,hetdat,hetvar);

fstar=delta'*rmat'*inv(rmat*vdel*rmat')*rmat*delta;
ftest=(bigt-(i+1)*q-p)*fstar/(bigt*i);

