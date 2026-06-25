% NAME            : Structural_Changes.m
% Description     :        
% Usage           : 
% Input           : 
% Input           : 
% Input           : 
%                                       
% Input           : sqrt(evar_m(JI)) - square root of right window spot variance estimate
% Input           : Test_Val - is value of test stastic
% Output          : beta_eff  = Jump Beta Efficient Estimator 
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
%
%
%==========================================================================
% COMMENTS
%==========================================================================

%% ========================================================================
%             Preliminaries
%  ========================================================================
% close all; clear all; clc;

yyy = csvread('real.csv');                                                        % read data
figure(1)
plot(yyy);
load('MarketData.mat');
yyy=MarketData(:,3);
T = length(yyy(1:800));                                                                % T = Sample size
bigt=T;                                                                       % set effective sample size
y=yyy(1:T,1);                                                                 % set up the data, y is the dependent variable z is 
                                                                                % the matrix of regressors (bigt,q) whose coefficients 
                                                                                % are allowed to change, x is a (bigt,p) matrix of 
                                                                                % regressors with coefficients fixed across regimes. 
                                                                                % Note: initialize x to [] if p=0.
% ------------
load('MSFTret.mat');                                                         
load('MarketData.mat');

A = MSFTret;                                                                    % Store returns in A and B
B = MarketData;
ExampleData = [B(:,3) A(:,3)];

%% ========================================================================
%             Initialize Seed
%  ========================================================================
rand('seed',19330407);                                                          %April 7, 1933
randn('seed',18650409);                                                         %April 9, 1865

X = ExampleData(:,1);                                                           % X is the asset
Y = ExampleData(:,2);                                                           % Y is the market
freq = 77;                                                                      % There are 77 5-minute returns in a normal trading day
kn = 15;                                                                        % Sub-sequence of time steps. It eliminates 15 intial and 15 end observation per each day? Because volatility is too high at the beginning anf end of the day 
T = floor(size(X,1)/freq);                                                      % # of days 
alpha_jump = 4.5;                                                               % Truncation for jump
alpha_cont = 3;

%% ========================================================================
%            Make space for Statistics
%  ========================================================================
D = zeros(5,12);                                                                % 12 outcome statistics, six years and full sample


%% ========================================================================
%            Compute Realized Variation and Bipower Variation and tod
%  ========================================================================
bv_X = zeros(T,1); 
rv_X = zeros(T,1);
bv_Y = zeros(T,1); 
rv_Y = zeros(T,1);
for t=1:T
    bv_X(t) = (freq/(freq-1))*(pi/2) * sum(abs(X((t-1)*freq+2:t*freq).*X((t-1)*freq+1:t*freq-1)));
    rv_X(t) = sum(X((t-1)*freq+1:t*freq).^2);
    bv_Y(t) = (freq/(freq-1))*(pi/2)*sum(abs(Y((t-1)*freq+2:t*freq).*Y((t-1)*freq+1:t*freq-1)));
    rv_Y(t) = sum(Y((t-1)*freq+1:t*freq).^2);
end
%figure();
%plot(rv_X(1:100));

figure()
plot(rv_Y)



% -----------
yyy = csvread('real.csv');                                                        % read data
figure(1)
plot(yyy);
load('MarketData.mat');
%yyy=jump;
T = length(yyy);                                                                % T = Sample size
bigt=T;                                                                       % set effective sample size
y=yyy(1:T,1);                                                                 % set up the data, y is the dependent variable z is 
                                                                                % the matrix of regressors (bigt,q) whose coefficients 
                                                                                % are allowed to change, x is a (bigt,p) matrix of 
                                                                                % regressors with coefficients fixed across regimes. 
                                                                                % Note: initialize x to [] if p=0.



z=ones(bigt,1);

yyy=newyy;
z=newx;
T = length(yyy);                                                                % T = Sample size
bigt=T;                                                                       % set effective sample size
y=yyy(1:T,1); 



x=[];

q=1;                                                                            % Number of regressors z potentially associated with structural chnage
p=0;                                                                            % Number of regressors x 
m=5; %default                                                                   % Maximum number of structural changes allowed. Defaul is m=5. But see below.
%  ========================================================================
%             Trimming Parameter (for supF Test)
%  ========================================================================

eps1=.05;                                                                       % value of the trimming (in percentage) for the construction 
                                                                                % and critical values of the supF type tests (used in the 
                                                                                % supF test, the Dmax, the supF(l+1|l) and the sequential 
                                                                                % procedure). If these tests are used, h below should be set 
                                                                                % at int(eps1*bigt). But if the tests are not required, 
                                                                                % estimation can be done with an arbitrary h. There are five 
                                                                                % options: eps1= .05, .10, .15, .20, or .25. 
                                                                                
                                                                                % For each option, 
                                                                                % the maximal value of m above is: 10 for eps=.05, 8 for 
                                                                                % eps1=.10, 5 for eps1=.15, 3 for eps1=.20, and 2 for eps1=.25.

h=round(eps1*bigt);                                                             % minimal length of a segment (h >= q). Note: if robust=1, h 
                                                                                % should be set at a larger value. 
%  ========================================================================
%             Options for Partial Change Model (p>0)
%  ========================================================================
fixb=0;                   % set to 1 if use fixed initial values for beta
betaini=0;                % if fixb=1, load the initial value of beta
maxi=20;                  % maximum number of iterations for the nonlinear procedure to 
                          % obtain global minimizers
printd=1;                 % set to 1 if want the output from the iterations to be printed
eps=0.0001;               % criterion for the convergence
% --------------------------------------------------------------------


%  ========================================================================
%            General Options 
%  ========================================================================
%        Option for Heterogeneity and Autocorrelation in Errors
%  ------------------------------------------------------------------------
robust=0;                                                                       % set to 1 if want to allow for heterogeneity and autocorrelation 
                                                                                % in the residuals, 0 otherwise. The method used is Andrews(1991) 
                                                                                % automatic bandwidth with AR(1) approximation and the quadratic 
                                                                                % kernel. Note: Do not set to 1 if lagged dependent variables are 
                                                                                % included as regressors.
% -------------------------------------------------------------------------
%             Option for pre-Whitening in the estimation of LR variance              
% -------------------------------------------------------------------------
prewhit=0;                                                                      % set to 1 if want to apply AR(1) prewhitening prior to estimating 
                                                                                % the long run covariance matrix.
                                                                                
% =========================================================================
%             Option for F tests   
% =========================================================================
% -------------------------------------------------------------------------
%             Regressors             
% -------------------------------------------------------------------------

hetdat=1;                                                                       % option for the construction of the F tests. Set to 1 if want to
                                                                                % allow different moment matrices of the regressors across segments. 
                                                                                % If hetdat=0, the same moment matrices are assumed for each segment 
                                                                                % and estimated from the ful sample. 
                                                                                % It is recommended to set hetdat=1 if p>0.                                                                                
% -------------------------------------------------------------------------
%             Variance Errors             
% -------------------------------------------------------------------------                                                                                
hetvar=1;                                                                       % option for the construction of the F tests.Set to 1 if want to allow 
                                                                                % for the variance of the residuals to be different across segments. 
                                                                                % If hetvar=0, the variance of the residuals is assumed constant 
                                                                                % across segments and constructed from the ful sample. This option 
                                                                                % is not available when robust=1.  
% =========================================================================
%             Option for CI in the break dates  
% =========================================================================                                                                                
hetomega=1;                                                                     % used in the construction of the confidence intervals for the break 
                                                                                % dates. If hetomega=0, the long run covariance matrix of zu is 
                                                                                % assumed identical across segments (the variance of the errors u 
                                                                                % if robust=0)
hetq=1;                                                                         % used in the construction of the confidence intervals for the break 
                                                                                % dates. If hetq=0, the moment matrix of the data is assumed identical 
                                                                                % across segments.
% =========================================================================
%             Global Minimizer Procedure  
% =========================================================================                                                                                
doglobal=1;                                                                     % set to 1 if want to cal the procedure to obtain global minimizers
dotest=1;                                                                       % set to 1 if want to construct the supF, UDmax and WDmax tests 
                                                                                % doglobal must be set to 1 to run this procedure.
dospflp1=1;                                                                     % set to 1 if want to construct the supF(l+1|l) tests where under
                                                                                % the null the l breaks are obtained using global minimizers. 
                                                                                % doglobal must be set to 1 to run this procedure.
% =========================================================================
%            Information Criterion  
% =========================================================================                                                                                
doorder=1;                                                                      % set to 1 if want to cal the procedure that selects the number of
                                                                                % breaks using information criteria. doglobal must be set to 1 to 
                                                                                % run this procedure.
% =========================================================================
%            Sequential Procedure 
% =========================================================================                                                                                
dosequa=1;                                                                      % set to 1 if want to estimate the breaks sequentialy and estimate 
                                                                                % the number of breaks using supF(l+1|l) test   
% =========================================================================
%            Option for Repartition Method of Bai (1995
% =========================================================================                                                                                
dorepart=1;                                                                     % set to 1 if want to modify the break dates obtained from the 
                                                                                % sequential method using the repartition method of Bai (1995),
                                                                                % Estimating breaks one at a time. This is needed for the confidence 
                                                                                % intervals obtained with estim below to be valid.
                                                                                
% =========================================================================
%            Selection Estimation Method for Number of Breaks
% =========================================================================                                                                                
estimbic=0;              % set to 1 if want to estimate the model with the number of breaks 
                         % selected by BIC.
estimlwz=0;              % set to 1 if want to estimate the model with the number of breaks  
                         % selected by LWZ
estimseq=1;              % set to 1 if want to estimate the model with the number of breaks
                         % selected using the sequential procedure
estimrep=1;              % set to 1 if want to estimate the model with the breaks selected
                         % using the repartition method
% -------------------------------------------------------------------------
%             Prespecified Number of Breaks            
% -------------------------------------------------------------------------                         
estimfix=0;                                                                     % set to 1 if want to estimate the model with a prespecified number
                                                                                % of breaks equal to fixn set below
fixn=0;

pbreak(bigt,y,z,q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,dotest,dospflp1,doorder,dosequa,dorepart,estimbic,estimlwz,estimseq,estimrep,estimfix,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar,fixn)
