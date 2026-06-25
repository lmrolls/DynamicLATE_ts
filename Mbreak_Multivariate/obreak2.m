
% This MATLAB version: January 17, 2022.
% This program is distributed freely for non-profit academic purposes only. for other uses, please contact Zhongjun Qu (zqu at uiuc.edu).
% A lot of effort has been put to construct this program and we would appreciate that
% you acknowledge using this code in your research and cite the relevant paper on which it is based: 
% Qu, Z. and P. Perron (2007): "Estimating and Testing for Structural Changes in Multivariate Regressions", Econometrica.
% Although a lot of efforts have been put in constructing the program, we cannot be
% held responsible for any consequences that could result from remaining errors.
% Copyright: Zhongjun Qu and Pierre Perron (2007)


clear
clc

% load the data

load data2.dat;
data = data2;

% Set the number of observations
global T             
T=100;     

% Set your T by n matrix of dependent variables         
                              
y=data(:,1:2);                           

% Set your T by q matrix of regressors  

z=data(:,3:5);                            

% Specify the matrix _S (here the first and second regressors are used in the 
% first equation and the first and third regressors are used in the second equation)   
global S_                                                                       
S_=[1 0 0 0                             
    0 1 0 0                             
    0 0 0 0  
    0 0 1 0 
    0 0 0 0 
    0 0 0 1];   

% Specify the number of regressors in the first (q1) and second (q2) equations.
                                    
q1=2;                                     
q2=2;                                     

% set the trimming parameter that specifes the minimal lenght of a segment as a proportion of the sample size.@

trm=0.2; 

% Set the trunction needed for constructing the test@

trun=floor(sqrt(T));

% set to getcon=1 if you want to construct the confidence intervals for the estimates of the break dates@                              

getcon=0; 

% Set hetq=1 if you want to allow the distribution of the regressors to change across
% regimes (this is used only when constructing confidence intervals for the estimates 
% of the break dates. For the construction of the test hetq=1 always.
global hetq 
hetq=1;                                  

% Set the bound for simulating the critical values, see the paper.

bigM=20;  

% Set the number of replications                              

rep=5000;            

% The current version does not allow restrictions so R is set to an identity matrix, R=eye(cols(_S)*(m+1)).@
global R
R=eye(cols(S_)*2);                   


% call the main procedure

mainop(y,z,cols(y),q1,q2,trm,T,hetq,trun,getcon,bigM,rep);




