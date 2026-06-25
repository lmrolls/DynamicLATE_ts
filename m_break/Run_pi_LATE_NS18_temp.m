% NAME            : Run_pi_LATE.m
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
restoredefaultpath
close all; clear all; clc;
dta = readtable('dataForMatlabAltEstTick.csv');
trim_date = datetime(datestr(datenum(dta.date_daily(401:end), 'ddmmyyyy')));
FOMC = dta.FOMCused(401:end);
DNY = dta.DNY2(401:end);
DNF = dta.DNF2(401:end);
DRF = dta.DRF2(401:end);
it = dta.path_intra_wide(401:end);
cont_date = trim_date(FOMC==0);
treat_date = trim_date(FOMC==1);
D = DNY;
y1 = DNF;
y2 = DRF;
T = size(DNY,1); % time dimension
Z = zeros(T,1); % instruments
T_C=size(cont_date,1);
T_P=size(treat_date,1);
Z(FOMC==1) = T/T_P;
Z(FOMC==0) = -T/T_C;
%Z = Z.*D;
D=D.^2;
Dc = D(FOMC==0);
Zc = Z(FOMC==0);
Dp = D(FOMC==1);
Zp = Z(FOMC==1);
desired_ci = 4;
D_for_Test=[Dc; Dp];
Z_for_Test=[Zc; Zp];
bigt=T;                                                                         % set effective sample size
y=D_for_Test;                                                                   % set up the data, y is the dependent variable z is 
                                                                                % the matrix of regressors (bigt,q) whose coefficients 
                                                                                % are allowed to change, x is a (bigt,p) matrix of 
                                                                                % regressors with coefficients fixed across regimes. 
                                                                                % Note: initialize x to [] if p=0.

z=Z_for_Test;

%x=ones(T,1);
x=[];

q=1;                                                                            % Number of regressors z potentially associated with structural chnage
p=0;                                                                            % Number of regressors x 
m=5; %default                                                                   % Maximum number of structural changes allowed. Defaul is m=5. But see below.

%  ========================================================================
%             Trimming Parameter (for supF Test)
%  ========================================================================
eps1 = .05;                                                                     % The asymptotic distribution depends on a trimming parameter
                                                                                % via the imposition of the minimal length h of a segment, namely
                                                                                % h = eps1*T     
                                                                                % value of the trimming (in percentage) for the construction 
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
robust=1;                                                                       % set to 1 if want to allow for heterogeneity and autocorrelation 
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
doglobal=1;                                                                     % set to 1 if want to call the procedure to obtain global minimizers
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
estimseq=0;              % set to 1 if want to estimate the model with the number of breaks
                         % selected using the sequential procedure
estimrep=0;              % set to 1 if want to estimate the model with the breaks selected
                         % using the repartition method
% -------------------------------------------------------------------------
%             Prespecified Number of Breaks            
% -------------------------------------------------------------------------                         
estimfix=0;                                                                     % set to 1 if want to estimate the model with a prespecified number
                                                                                % of breaks equal to fixn set below
fixn=0;

Procedure_pi_LATE(bigt,y,z,q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,dotest,dospflp1,doorder,dosequa,dorepart,estimbic,estimlwz,estimseq,estimrep,estimfix,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar,fixn)

T=bigt;
q=size(z,2);
p=size(x,2);
if isempty(x)=1
[glob,datevec,bigvec]=Procedure_piL_Full_Sample(bigt,y(1:T_C),z(1:T_C),q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar);
else
    [glob,datevec,bigvec]=Procedure_piL_Full_Sample(bigt,y(1:T_C),z(1:T_C),q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,fixb,x(1:T_C),p,eps,maxi,betaini,printd,hetdat,hetvar);
end
%%=========================================================================
%                               ST Dates
%==========================================================================
%vector_piL=[0.4:0.1:1];
vector_piL=0.6;
Lenght_piL=length(vector_piL);
%max_breaks=5;
max_breaks=1;
step_pi=0.05;
for piL=1:Lenght_piL
num_pis=vector_piL(piL):step_pi:1;
index_pis=length(num_pis);
F_ST=zeros(index_pis,max_breaks);
for m=1:max_breaks
%for m=5:5
for pi=1:index_pis
    S_index=zeros(T,1);
    Tpi=floor(num_pis(pi)*T);
    if vector_piL(piL)==1
        S_m=[1,T];
    else
    %S_m=S_Dates_m(y, z, x, datevec, T, num_pis(pi), m);
    S_m=S_Dates_m_Heter(y, z, x, datevec, T, num_pis(pi), m);
    end
    rows = size(S_m,1);
        for irow=1:rows
            S_index(S_m(irow,1):S_m(irow,2))=1;
        end

y_ST=y(S_index==1);
z_ST=z(S_index==1);
if isempty(x)==0
    x_ST=x(S_index==1);
    Mx_ST=x_ST*(x_ST'*x_ST)^-1*x_ST';
else
    x=zeros(T,1);
    x_ST=x(S_index==1);
    Mx_ST=eye(Tpi);
    x=[];
end
y_tilde=Mx_ST*y_ST;
z_tilde=Mx_ST*z_ST;
F_ST(pi,m)=F_ST_Stat(y_tilde, z_tilde, Tpi, q, p);
end
end
[Val_m, ind_m]=max(F_ST);
max_pi_for_m=num_pis(ind_m);
[Val_m_pi, ind_m_pi]=max(Val_m);
Sup_m=ind_m_pi;
Sup_pi=num_pis(ind_m(Sup_m));
[Sup_F, maxIndex]=max(F_ST, [], 'all');
ind_all=find(max(F_ST));
disp(['pi_L:', num2str(vector_piL(piL))]);
disp(['Maximum value:', num2str(Sup_F)]); 
disp(['Maximum m:', ' ', num2str(Sup_m)]);
disp(['Maximum pi:', ' ', num2str(Sup_pi)]);

Sup_F_piL(piL,1)=Sup_F;
Sup_m_piL(piL,1)=Sup_m;
Sup_pi_piL(piL,1)=Sup_pi;

end

Sup_S_m=S_Dates_m(y, z, x, datevec, T, Sup_pi, Sup_m);
Index_S=zeros(T,1);
rows = size(Sup_S_m,1);
        for irow=1:rows
            Index_S(Sup_S_m(irow,1):Sup_S_m(irow,2))=1;
        end
y_S=y(Index_S==1);
z_S=z(Index_S==1);
if isempty(x)==0
    x_S=x(Index_S==1);
else
    x_S=x;
end
%%=========================================================================
%                   Estimation in pi Late sample
%==========================================================================
[estimator, t_stat, Fstat, Sup_S_m]=Estimation_pi_Late(y_S, z_S, x_S, T, Sup_S_m, Sup_pi);