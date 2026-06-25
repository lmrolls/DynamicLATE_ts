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
%Y_Nominal = DNF;
%Y_Real = DRF;
T = size(DNY,1); % time dimension
Z = zeros(T,1); % instruments
T_C=size(cont_date,1);
T_P=size(treat_date,1);
Z(FOMC==1) = T/T_P;
Z(FOMC==0) = -T/T_C;
Z=Z.*D;
Dc = D(FOMC==0);
Dp = D(FOMC==1);
D_for_Test=[Dc-mean(Dc,1); Dp-mean(Dp,1)];
Zc = -T/T_C*(Dc-mean(Dc,1));
Zp = T/T_P*(Dp-mean(Dp,1));
%Z = Z.*D;
%D=D.^2;
%Dc = D(FOMC==0);
%Zc = Z(FOMC==0);
%Zp = Z(FOMC==1);
%desired_ci = 4;
%D_for_Test=[Dc-mean(Dc); Dp-mean(Dp)];
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

h=round(eps1*T_C);                                                             % minimal length of a segment (h >= q). Note: if robust=1, h 
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
if isempty(x)==1
[glob,datevec,bigvec]=Procedure_piL_Full_Sample(T_C,y(1:T_C),z(1:T_C),q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar);
else
    [glob,datevec,bigvec]=Procedure_piL_Full_Sample(T_C,y(1:T_C),z(1:T_C),q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,fixb,x(1:T_C),p,eps,maxi,betaini,printd,hetdat,hetvar);
end
clear m
%%=========================================================================
%                               ST Dates
%==========================================================================
Hac_option=0;                                                               % Hac_option=0 for heteroskedasticity-only. Hac_option=1 for HAC.   
%vector_piL=[0.6:0.1:1];
vector_piL=0.8;
Lenght_piL=length(vector_piL);
max_breaks=5;
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
    S_m=S_Dates_m_Heter(y, z, x, datevec, T, num_pis(pi), m, T_C);
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
F_ST(pi,m)=F_ST_Stat(y_tilde, z_tilde, Tpi, q, p, Hac_option);
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

%Sup_S_m=S_Dates_m(y, z, x, datevec, T, Sup_pi, Sup_m);                     %Not fixing the treatment sample 
Sup_S_m=S_Dates_m_Heter(y, z, x, datevec, T, Sup_pi, Sup_m, T_C);           %Fixing the treatment sample. Whole P sample selected
%disp(['pi_L:', num2str(vector_piL(piL))]);
disp(['Srongly-Identified S:']);
Sup_S_m

end

[value, index_Sup_piL]=max(Sup_F_piL);
%Sup_S_m=S_Dates_m(y, z, x, datevec, T, Sup_pi, Sup_m);                     %Not fixing the treatment sample 
Sup_S=S_Dates_m_Heter(y, z, x, datevec, T, Sup_pi_piL(index_Sup_piL,1), Sup_m_piL(index_Sup_piL,1), T_C);           %Fixing the treatment sample. Whole P sample selected
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
y_S=y(Index_S==1);
z_S=z(Index_S==1);
if isempty(x)==0
    x_S=x(Index_S==1);
else
    x_S=x;
end
pi_hat=Sup_pi_piL(index_Sup_piL,1);

DTA = readtable('results/NY.csv');
windows = readtable('results/windows.csv');
load('yearscount/time.mat')
D = DTA.Yields; 
FOMC = DTA.FOMC;
D_C = D(FOMC==0);
D_P = D(FOMC==1);
index_larger = find(Sup_S(1:end,:) > T_C, 1, 'first');
for i=1:size(Sup_S(1:index_larger-1,:),1)
    window_dates(2*i-1:2*i,1)=Sup_S(i,:)';
end
% Control Sample
figure()
subplot(2,1,1)
plot(timeC,D_C,'LineWidth',1.2) %series
for j = 1:(length(window_dates)/2)
rectangle('Position',[timeC(window_dates(2*j-1)) min(D_C) timeC(window_dates(2*j)) - timeC(window_dates(2*j-1)) max(D_C)-min(D_C)],'LineStyle','--','LineWidth',1.2,'EdgeColor','r') %window 1
text(timeC(round(window_dates(2*j)-(window_dates(2*j)-window_dates(2*j-1))/2)),max(D_C)+0.055,['${S_{', num2str(j), '}}$' ], 'FontWeight','bold', 'Interpreter','latex')
end
axis([timeC(1) timeC(end) -0.35 0.35])
lgd = legend;
lgd.FontSize = 8;
lgd.Title.String = ['$\hat{\pi}=', num2str(pi_hat),'$']; 
lgd.String='$\Delta_{C}$';
lgd.Position=[0.785595237414043,0.82976753241888,0.119642857824053,0.095238097054618];
lgd.Interpreter='Latex';
title('Control Sample')

% Treatment Sample
for i=1:size(Sup_S(index_larger:end,:),1)
    window_dates_P(2*i-1:2*i,1)=Sup_S(index_larger+i-1,:)'-T_C;
end
subplot(2,1,2)
plot(timeP,D_P,'LineWidth',1.2) %series
for j = 1:(length(window_dates_P)/2)
rectangle('Position',[timeP(window_dates_P(2*j-1)) min(D_P) timeP(window_dates_P(2*j)) - timeP(window_dates_P(2*j-1)) max(D_P)-min(D_P)],'LineStyle','--','LineWidth',1.2,'EdgeColor','r') %window 1
text(timeP(round(window_dates_P(2*j)-(window_dates_P(2*j)-window_dates_P(2*j-1))/2)),max(D_P)+0.055,['${S_{', num2str(j), '}}$' ], 'FontWeight','bold', 'Interpreter','latex')
end
axis([timeP(1) timeP(end) -0.35 0.35])
lgd = legend;
lgd.FontSize = 8;
lgd.Title.String = ['$\hat{\pi}=', num2str(pi_hat),'$']; 
lgd.String='$\Delta_{P}$';
lgd.Position=[0.785595237414044,0.355481818116947,0.119642857824053,0.095238097054618];
lgd.Interpreter='Latex';
title('Treatment Sample')
%saveas(gcf,'results\figures\Figure2.eps','epsc')


%%=========================================================================
%                   Estimation of tetha in pi Late sample
%==========================================================================
[theta_hat, t_stat, Fstat, Sup_S]=Estimation_pi_Late(y_S, z_S, x_S, T, Sup_S, pi_hat, Hac_option);
%[estimator, t_stat1, Fstat, Sup_S, SE_robust, CI]=Estimation_pi_Late(Y, D, x, T, Sup_S, pi, HAC)

% 2SLS Estimation
Y_Nominal = DNF;
Y_Real = DRF;
Y_Nominal_C = Y_Nominal(FOMC==0);
Y_Nominal_P = Y_Nominal(FOMC==1);
Y_Real_C = Y_Real(FOMC==0);
Y_Real_P = Y_Real(FOMC==1);
Y_Nominal_parti=[Y_Nominal_C;Y_Nominal_P];
Y_Real_parti=[Y_Real_C;Y_Real_P];
Y_Nom_S=Y_Nominal_parti(Index_S==1);
Y_Real_S=Y_Real_parti(Index_S==1);
D_S=y_S;
Z_S=z_S;
P_Z_S=Z_S*inv(Z_S'*Z_S)*Z_S';
D_hat_S=P_Z_S*D_S;
[beta_hat_S, t_stat1_S, Fstat, Sup_S, SE_robust_S, CI_S]=Estimation_pi_Late(Y_Real_S, D_hat_S, x_S, T, Sup_S, pi_hat, Hac_option);
[beta_hat_S_nom, t_stat1_S_nom, Fstat_nom, Sup_S, SE_robust_S_nom, CI_S_nom]=Estimation_pi_Late(Y_Nom_S, D_hat_S, x_S, T, Sup_S, pi_hat, Hac_option);


%%=========================================================================
%                   Full-Sample Estimation 
%==========================================================================
%FOMC = dta.FOMCused(401:end);
%DNY = dta.DNY2(401:end);
%DNF = dta.DNF2(401:end);
%DRF = dta.DRF2(401:end);
%it = dta.path_intra_wide(401:end);
%cont_date = trim_date(FOMC==0);
%treat_date = trim_date(FOMC==1);
% D = DNY;
% Y_Nominal = DNF;
% Y_Real = DRF;
% T = size(DNY,1); % time dimension
% Z = zeros(T,1); % instruments
% T_C=size(cont_date,1);
% T_P=size(treat_date,1);
% Z(FOMC==1) = T/T_P;
% Z(FOMC==0) = -T/T_C;
% Dc = D(FOMC==0);
% Dp = D(FOMC==1);
% D_for_Test=[Dc-mean(Dc,1); Dp-mean(Dp,1)];
% Zc = -T/T_C*(Dc-mean(Dc,1));
% Zp = T/T_P*(Dp-mean(Dp,1));
% Z_for_Test=[Zc; Zp];
% y=D_for_Test;                                                                   
% z=Z_for_Test;


D_parti=y;
Z_parti=z;
P_Z_parti=Z_parti*inv(Z_parti'*Z_parti)*Z_parti';
D_hat_parti=P_Z_parti*D_parti;

P_Z=Z*inv(Z'*Z)*Z';
D_hat=P_Z*D;

%% Real Forward Rates
[beta_hat_full, t_stat1_beta_full, Fstat_beta_full, RSS,  SE_robust_beta_full, CI_beta_full] = OLS_Robust_full(Y_Real_parti,D_hat_parti,x, 0, Hac_option);
%[beta_hat_full, t_stat1_beta_full, Fstat_beta_full, RSS,  SE_robust_beta_full, CI_beta_full] = OLS_Robust_full(Y_Real,D_hat,x, 0, Hac_option);

%NS Estimator
Cov_P_YD=cov(Y_Real(FOMC==1), D(FOMC==1)) ;
Cov_C_YD=cov(Y_Real(FOMC==0), D(FOMC==0));
Cov_P_YD=Cov_P_YD(1,2);
Cov_C_YD=Cov_C_YD(1,2);
beta_hat_NS=(Cov_P_YD - Cov_C_YD) / (var(D(FOMC==1)) - var(D(FOMC==0)));

%% Nominal Forward Rates
[beta_hat_full_nom, t_stat1_beta_full_nom, Fstat_beta_full_nom, RSS,  SE_robust_beta_full_nom, CI_beta_full_nom] = OLS_Robust_full(Y_Nominal_parti,D_hat_parti,x, 0, Hac_option);
%[beta_hat_full, t_stat1_beta_full, Fstat_beta_full, RSS,  SE_robust_beta_full, CI_beta_full] = OLS_Robust_full(Y_Real,D_hat,x, 0, Hac_option);

%NS Estimator, Nominal
Cov_P_YD_nom=cov(Y_Nominal(FOMC==1), D(FOMC==1)) ;
Cov_C_YD_nom=cov(Y_Nominal(FOMC==0), D(FOMC==0));
Cov_P_YD_nom=Cov_P_YD_nom(1,2);
Cov_C_YD_nom=Cov_C_YD_nom(1,2);
beta_hat_NS_nom=(Cov_P_YD_nom - Cov_C_YD_nom) / (var(D(FOMC==1)) - var(D(FOMC==0)));




















