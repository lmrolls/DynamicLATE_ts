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
% close all; clear all; clc;

function [S_m]=S_Dates_m_M2_Heter(Y, D, Z, X, Mx, beta_H0, datevec, T, pi,m, T_C)
Tpi=floor(pi*T);
vec_M2_ST=zeros(T,1);
q=size(Z,2);
p=size(X,2);
%%=========================================================================
%                               m=1 
%==========================================================================
if m==1
    Break_Date=[1, datevec(1,1);
                datevec(1,1)+1, T_C];
elseif m==2
    Break_Date=[1, datevec(1,2);
                datevec(1,2)+1, datevec(2,2);
                datevec(2,2)+1, T_C];
elseif m==3
    Break_Date=[1, datevec(1,3);
                datevec(1,3)+1, datevec(2,3);
                datevec(2,3)+1, datevec(3,3);
                datevec(3,3)+1, T_C];
elseif m==4
    Break_Date=[1, datevec(1,4);
                datevec(1,4)+1, datevec(2,4);
                datevec(2,4)+1, datevec(3,4);
                datevec(3,4)+1, datevec(4,4);
                datevec(4,4)+1, T_C];
elseif m==5
    Break_Date=[1, datevec(1,5);
                datevec(1,5)+1, datevec(2,5);
                datevec(2,5)+1, datevec(3,5);
                datevec(3,5)+1, datevec(4,5);
                datevec(4,5)+1, datevec(5,5);
                datevec(5,5)+1, T_C];
end
for i=1:m+1
    CT_Z_S=Z(Break_Date(i,1):Break_Date(i,2),:);
    T_S=size(CT_Z_S,1);
    if isempty(X)==0
        X_S=X(Break_Date(i,1):Break_Date(i,2),:);
        Mx_S=eye(T_S)-X_S*inv(X_S'*X_S)*X_S';
    else
        Mx_S=eye(T_S);
    end
    Zbar_S=Mx_S*CT_Z_S;
    %Zbar_S=Mx*CT_Z_S;
    y_S=[Y(Break_Date(i,1):Break_Date(i,2),1), D(Break_Date(i,1):Break_Date(i,2),1)];
    a=[beta_H0;1];
    b=[1; -beta_H0];
    % Statistics N1 and N2
    %v*b.*Zbar;
     if isempty(X)==0
         v_hat_S=y_S-Zbar_S*inv(Zbar_S'*Zbar_S)*Zbar_S'*y_S- (eye(size(Mx_S,1))-Mx_S)*y_S;
     else
         v_hat_S=y_S-Zbar_S*inv(Zbar_S'*Zbar_S)*Zbar_S'*y_S;
     end
Sigma_v_hat=v_hat_S'*v_hat_S/(T_S-q-p);
Lag=floor(T_S^(1/3));
Sigma_N1_hat =  NeweyWest(v_hat_S*b, Zbar_S,Lag,0);
Sigma_N2star_hat =  NeweyWest(v_hat_S*inv(Sigma_v_hat)*a, Zbar_S,Lag,0);
Sigma_N1N2_hat=NeweyWest_N1N2((v_hat_S*b.*Zbar_S), (v_hat_S*inv(Sigma_v_hat)*a.*Zbar_S),Lag,0);
Sigma_N2_hat=Sigma_N2star_hat-Sigma_N1N2_hat*inv(Sigma_N1_hat) *Sigma_N1N2_hat';
N1_ST=Sigma_N1_hat^(-1/2)*T_S^(-1/2)*Zbar_S'*y_S*b;
N2_ST=Sigma_N2_hat^(-1/2)*(T_S^(-1/2)*Zbar_S'*y_S*inv(Sigma_v_hat)*a-Sigma_N1N2_hat*Sigma_N1_hat^(-1/2)*N1_ST);
M2_ST(i,1)=N2_ST'*N2_ST;
Table_t_stat_dates(i,1:3)=[M2_ST(i,1), Break_Date(i,1), Break_Date(i,2)];
end
Table_t_stat_dates=[Table_t_stat_dates; inf, T_C+1, T];
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)+1];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>=Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
if ind_smaller<1
    Missing=Tpi;
    Start_date=Table_t_stat_dates(ind_larger,2)-1;
    End_date=Table_t_stat_dates(ind_larger,3);
else
    Missing=Tpi-A(ind_smaller,1);
    Start_date=Table_t_stat_dates(ind_larger,2)-1;
    End_date=Table_t_stat_dates(ind_larger,3);
end
Length_reg=End_date-Start_date;
if Length_reg==Missing 
    if m==1
        pi_Dates=[1,T];
    else
        if ind_smaller<1
            if Start_date==0
                Start_date=1;
            end
            pi_Dates=[Start_date,End_date-1];
        else     
            if Start_date==0
                Start_date=1;
            end
            pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+1,End_date];
            pi_Dates=sortrows(pi_Dates, 1,'ascend');
        end
    end
else
for i=1:Length_reg-Missing
    CT_Z_S=Z(Start_date+i:Start_date + i +Missing,:);
    T_S=size(CT_Z_S,1);
    if isempty(X)==0
        X_S=X(Start_date+i:Start_date + i +Missing,:);
        Mx_S=eye(T_S)-X_S*inv(X_S'*X_S)*X_S';
    else
        Mx_S=eye(T_S);
    end
    Zbar_S=Mx_S*CT_Z_S;
    y_S=[Y(Start_date+i:Start_date + i +Missing,1), D(Start_date+i:Start_date + i +Missing,1)];
    a=[beta_H0;1];
    b=[1; -beta_H0];
    % Statistics N1 and N2
    %v*b.*Zbar;
     if isempty(X)==0
         v_hat_S=y_S-Zbar_S*inv(Zbar_S'*Zbar_S)*Zbar_S'*y_S- (eye(size(Mx_S,1))-Mx_S)*y_S;
     else
         v_hat_S=y_S-Zbar_S*inv(Zbar_S'*Zbar_S)*Zbar_S'*y_S;
     end
Sigma_v_hat=v_hat_S'*v_hat_S/(T_S-q-p);
Lag=floor(T_S^(1/3));
Sigma_N1_hat =  NeweyWest(v_hat_S*b, Zbar_S,Lag,0);
Sigma_N2star_hat =  NeweyWest(v_hat_S*inv(Sigma_v_hat)*a, Zbar_S,Lag,0);
Sigma_N1N2_hat=NeweyWest_N1N2((v_hat_S*b.*Zbar_S), (v_hat_S*inv(Sigma_v_hat)*a.*Zbar_S),Lag,0);
Sigma_N2_hat=Sigma_N2star_hat-Sigma_N1N2_hat*inv(Sigma_N1_hat) *Sigma_N1N2_hat';
N1_ST=Sigma_N1_hat^(-1/2)*T_S^(-1/2)*Zbar_S'*y_S*b;
N2_ST=Sigma_N2_hat^(-1/2)*(T_S^(-1/2)*Zbar_S'*y_S*inv(Sigma_v_hat)*a-Sigma_N1N2_hat*Sigma_N1_hat^(-1/2)*N1_ST);
vec_M2_ST(i,1)=N2_ST'*N2_ST;
end
[~, ind]=max(abs(vec_M2_ST));
    pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
    pi_Dates=sortrows(pi_Dates, 1,'ascend');
end
pi_Dates_m1=pi_Dates;
S_m1=pi_Dates_m1;
S_m=S_m1;
vec_M2_ST=zeros(T,1);
end