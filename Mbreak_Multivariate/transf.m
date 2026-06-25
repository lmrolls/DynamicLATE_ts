% @*************************************************************************@

function [maty,matz]=transf(y,z)

% @ 
% ***Procedure that generates the matrix of regressors (matz) and that of 
%    the dependent variable (maty)
% ***Input: y, dependent variable as a matrix, with the jth row containing 
%           observations at time j
%           z: regressors as a matrix, with the jth row containing 
%           observations at time j
% ***Output: maty, dependent variable as a row vector, with the (i-1)*n+jth 
%            row containing ith observation for equation j
%            matz: regressors as a matrix
% ***Global: _S, matrix of constants selecting variables appearing in each 
%           regression.
% @
global S_
% local maty,matz,i,_ID;
ID_=eye(cols(y));
maty=y(1,:)';
matz=kron(ID_,z(1,:))*S_;
i=2;
while i<=rows(y)
    maty=[maty; y(i,:)'];
    matz=[matz;kron(ID_,z(i,:))*S_];
    i=i+1;
end
% retp(maty,matz);
% endp;