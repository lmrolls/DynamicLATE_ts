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
 close all; clear all; clc;

yyy = csvread('real.csv');                                                        % read data
T = length(yyy);                                                                % T = Sample size                                                                        % set effective sample size
D=yyy(1:T,1);                                                                   % set up the data, y is the dependent variable z is 
                                                                                % the matrix of regressors (bigt,q) whose coefficients 
                                                                                % are allowed to change, x is a (bigt,p) matrix of 
                                                                                % regressors with coefficients fixed across regimes. 
                                                                                % Note: initialize x to [] if p=0.

Z=ones(T,1);

X=[];

q=1;                                                                            % Number of regressors z potentially associated with structural chnage
p=0;                                                                            % Number of regressors x 
m_breaks=5; %default                                                                   % Maximum number of structural changes allowed. Defaul is m=5. But see below.
Y=D+randn;
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

h=round(eps1*T);                                                             % minimal length of a segment (h >= q). Note: if robust=1, h 
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

%Procedure_pi_LATE(T,D,Z,q,m_breaks,h,eps1,robust,prewhit,hetomega,hetq,doglobal,dotest,dospflp1,doorder,dosequa,dorepart,estimbic,estimlwz,estimseq,estimrep,estimfix,fixb,X,p,eps,maxi,betaini,printd,hetdat,hetvar,fixn)

q=size(Z,2);
p=size(X,2);

[glob,datevec,bigvec]=Procedure_piL_Full_Sample(T,D,Z,q,m_breaks,h,eps1,robust,prewhit,hetomega,hetq,doglobal,fixb,X,p,eps,maxi,betaini,printd,hetdat,hetvar);

%%=========================================================================
%                               ST Dates
%==========================================================================
vector_piL=[0.6:0.1:1];
Lenght_piL=length(vector_piL);
for piL=1:Lenght_piL
max_breaks=5;
step_pi=0.05;
num_pis=vector_piL(piL):step_pi:1;
index_pis=length(num_pis);
F_ST=zeros(index_pis,max_breaks);
for m=1:max_breaks   
for pi=1:index_pis
    S_index=zeros(T,1);
    Tpi=floor(num_pis(pi)*T);
    S_m=S_Dates_m(D, Z, X, datevec, T, num_pis(pi), m);
    rows = size(S_m,1);
        for irow=1:rows
            S_index(S_m(irow,1):S_m(irow,2))=1;
        end

D_ST=D(S_index==1);
Z_ST=Z(S_index==1);
if isempty(X)==0
    X_ST=X(S_index==1);
    Mx_ST=X_ST*(X_ST'*X_ST)^-1*X_ST';
else
    X=zeros(T,1);
    X_ST=X(S_index==1);
    Mx_ST=eye(Tpi);
    X=[];
end
D_tilde=Mx_ST*D_ST;
Z_tilde=Mx_ST*Z_ST;
F_ST(pi,m)=F_ST_Stat(D_tilde, Z_tilde, Tpi, q, p);
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

Sup_S_m=S_Dates_m(D, Z, X, datevec, T, Sup_pi, Sup_m);                     %Not fixing the treatment sample 
%Sup_S_m=S_Dates_m_Heter(y, z, x, datevec, T, Sup_pi, Sup_m, T_C);           %Fixing the treatment sample. Whole P sample selected
%disp(['pi_L:', num2str(vector_piL(piL))]);
disp(['Srongly-Identified S:']);
Sup_S_m

end

[value, index_Sup_piL]=max(Sup_F_piL);
Sup_S=S_Dates_m(D, Z, X, datevec, T, Sup_pi_piL(index_Sup_piL,1), Sup_m_piL(index_Sup_piL,1));                     %Not fixing the treatment sample 
%Sup_S=S_Dates_m_Heter(y, z, x, datevec, T, Sup_pi_piL(index_Sup_piL,1), Sup_m_piL(index_Sup_piL,1), T_C);           %Fixing the treatment sample. Whole P sample selected
disp(['pi_L:', num2str(vector_piL(index_Sup_piL))]);
disp(['Maximum value:', num2str(Sup_F_piL(index_Sup_piL))]); 
disp(['Maximum m:', ' ', num2str(Sup_m_piL(index_Sup_piL,1))]);
disp(['Srongly-Identified S:']);
Sup_S

Index_S=zeros(T,1);
rows = size(Sup_S,1);
        for irow=1:rows
            Index_S(Sup_S(irow,1):Sup_S(irow,2))=1;
        end
D_S=D(Index_S==1);
Z_S=Z(Index_S==1);
if isempty(X)==0
    X_S=X(Index_S==1);
else
    X_S=X;
end
pi_hat=Sup_pi_piL(index_Sup_piL,1);
m_hat=num2str(Sup_m_piL(index_Sup_piL,1));
Tpi_hat=T*pi_hat;
%%=========================================================================
%                   Estimation in pi Late sample
%==========================================================================
Hac_option=0;                                                                   % Hac_option=0 for heteroskedasticity-only. Hac_option=1 for HAC.   
[theta_hat, t_stat_theta, Fstat, Sup_S]=Estimation_pi_Late(D_S, Z_S, X_S, T, Sup_S, pi_hat, Hac_option);

for i=1:size(Sup_S(1:end,:),1)
    window_dates(2*i-1:2*i,1)=Sup_S(i,:)';
end

figure()
plot(D,'LineWidth',1.2) %series
for j = 1:(length(window_dates)/2)
rectangle('Position',[(window_dates(2*j-1)) min(D) (window_dates(2*j) - window_dates(2*j-1)) max(D)-min(D)],'LineStyle','--','LineWidth',1.2,'EdgeColor','r') %window 1
text(round(window_dates(2*j)-(window_dates(2*j)-window_dates(2*j-1))/2),max(D)+0.055,['${S_{', num2str(j), '}}$' ], 'FontWeight','bold', 'Interpreter','latex')
end
axis([1 T min(D)-0.5 max(D)+0.5])
lgd = legend;
lgd.FontSize = 8;
lgd.Title.String = ['$\hat{\pi}=', num2str(pi_hat),'$']; 
lgd.String='$\Delta_{C}$';
lgd.Position=[0.785595237414043,0.82976753241888,0.119642857824053,0.095238097054618];
lgd.Interpreter='Latex';
title('Full Sample')


S_OLS_hat=Sup_S;
Index_S_OLS=Index_S;
%%=========================================================================
%                   Estimation of beta(S_OLS_hat)
%==========================================================================
Hac_option=0;
Y_S_OLS_hat=Y(Index_S_OLS==1);
D_S_OLS_hat=D(Index_S_OLS==1);
if isempty(X)==1
    X_S_OLS_hat=[];
else
    X_S_OLS_hat=X(Index_S_OLS==1);
end

[beta_hat_S_OLS, t_stat1_beta_S_OLS, Fstat_beta_S_OLS, RSS,  SE_robust_beta_S_OLS, CI_beta_S_OLS] = OLS_Robust_full(Y_S_OLS_hat,D_S_OLS_hat,X_S_OLS_hat, 0, Hac_option);




