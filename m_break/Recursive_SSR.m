function vec_Recursive_SSR = Recursive_SSR(start,y,z,h,last)

% NAME            : Recursive_SSR.m
% Description     : Computes recursive residuals from a data set that starts
%                   at date "start" and ends at date "last". It returns a vector of sum of
%                   squared residuals (SSR) of length last-start+1 (stored for convenience in
%                   a vector of length T).       
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

% ========================================================================
%              Initialize the Vectors
%  ========================================================================

vec_Recursive_SSR=zeros(last,1);                                                % vec_Recursive_SSR is the output vector

%% ========================================================================
%              Initialize the Recursion
%  ========================================================================
inv1 = inv(z(start:start+h-1,:)'*z(start:start+h-1,:));                         % 1st segment goes from start to start+h-1
                                                                                % z is  Txq

delta1=inv1*(z(start:start+h-1,:)'*y(start:start+h-1,1));                       % OLS estimator of 1st segment
res=y(start:start+h-1,1)-z(start:start+h-1,:)*delta1;                           % OLS residuals of 1st segment
vec_Recursive_SSR(start+h-1,1)=res'*res;                                        % ssr of recursive residuals of 1st segment

%% ========================================================================
%              Main Loop for Recursive Resdusals and Updating
%  ========================================================================
r = start + h;                                                                  % 2nd segment start at start + h
                                                                                % Note here we use while instead of standard loop
while r <= last
    v = y(r,1)-z(r,:)*delta1;
    invz = inv1*z(r,:)';                                                        % invz is (start+h)x1
    f = 1+z(r,:)*invz;                                                          % f is a scalar    
    delta2 = delta1+(invz*v)/f;                                                 % Update OLS estiamtor
    inv2 = inv1-(invz*invz')/f;                                                 % Update invz1 : Z'Z    
    inv1 = inv2;                                                                  
    delta1 = delta2;
    vec_Recursive_SSR(r,1) = vec_Recursive_SSR(r-1,1) + v*v/f;                  % Update Recursive residuals
    r=r+1;                                                                      % update index    
end

    