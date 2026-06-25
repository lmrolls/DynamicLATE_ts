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

function [S_m1, S_m2, S_m3, S_m4, S_m5]=S_Dates(y, z, x, datevec, T, pi)
Tpi=floor(pi*T);
vec_RSS=zeros(T,1);
%%=========================================================================
%                               m=1 
%==========================================================================
Break_Date=datevec(1,1);
if isempty(x)==1
    [estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Date),z(1:Break_Date,:),[], 0);
    [estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Date+1:T),z(Break_Date+1:T,:),[], 0);
else 
    [estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Date),z(1:Break_Date,:),x(1:Break_Date,:), 0);
    [estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Date+1:T),z(Break_Date+1:T,:),x(Break_Date+1:T,:), 0);
end

if size(z,2)==1
    Table_t_stat_dates=abs([t_stat1, 1, Break_Date(1,1);
        t_stat2, Break_Date(1,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
else
    Table_t_stat_dates=abs([Fstat1, 1, Break_Dates(1,1);
        Fstat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        Fstat3, Break_Dates(2,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
end


Missing=Tpi-A(ind_smaller,1)-1;
Start_date=Table_t_stat_dates(ind_larger,2);
End_date=Table_t_stat_dates(ind_larger,3);
Length_reg=End_date-Start_date;
for i=1:Length_reg-Missing
    if isempty(x)==1
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),[], 0);
    else
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),x(Start_date+i:Start_date + i +Missing,:), 0);
    end
end
    [value, ind]=min(vec_RSS);
    pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
    pi_Dates=sortrows(pi_Dates, 1,'ascend');
    
pi_Dates_m1=pi_Dates;
S_m1=pi_Dates_m1;
vec_RSS=zeros(T,1);
%%=========================================================================
%                               m=2 
%==========================================================================
Break_Dates=datevec(1:2,2);
if isempty(x)==1
[estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),[], 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),[], 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:T),z(Break_Dates(2,1)+1:T,:),[], 0);
else
    [estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),x(1:Break_Dates(1,1),:), 0);
    [estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),x(Break_Dates(1,1)+1:Break_Dates(2,1),:), 0);
    [estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:T),z(Break_Dates(2,1)+1:T,:),x(Break_Dates(2,1)+1:T,:), 0);
end


if size(z,2)==1
    Table_t_stat_dates=abs([t_stat1, 1, Break_Dates(1,1);
        t_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        t_stat3, Break_Dates(2,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
else
    Table_t_stat_dates=abs([Fstat1, 1, Break_Dates(1,1);
        Fstat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        Fstat3, Break_Dates(2,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
end


Missing=Tpi-A(ind_smaller,1)-1;
Start_date=Table_t_stat_dates(ind_larger,2);
End_date=Table_t_stat_dates(ind_larger,3);
Length_reg=End_date-Start_date;
for i=1:Length_reg-Missing
    if isempty(x)==1
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),[], 0);
    else
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),x(Start_date+i:Start_date + i +Missing,:), 0);
    end
end
    [value, ind]=min(vec_RSS);
    pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
    pi_Dates=sortrows(pi_Dates, 1,'ascend');
    pi_Dates_m2=pi_Dates;
    S_m2=pi_Dates;
vec_RSS=zeros(T,1);


%%=========================================================================
%                               m=3 
%==========================================================================
Break_Dates=datevec(1:3,3);
if isempty(x)==1
[estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),[], 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),[], 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:Break_Dates(3,1)),z(Break_Dates(2,1)+1:Break_Dates(3,1),:),[], 0);
[estimator4, t_stat4, Fstat4] = OLS_Robust(y(Break_Dates(3,1)+1:T),z(Break_Dates(3,1)+1:T,:),[], 0);
else
    [estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),x(1:Break_Dates(1,1),:), 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),x(Break_Dates(1,1)+1:Break_Dates(2,1),:), 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:Break_Dates(3,1)),z(Break_Dates(2,1)+1:Break_Dates(3,1),:),x(Break_Dates(2,1)+1:Break_Dates(3,1),:), 0);
[estimator4, t_stat4, Fstat4] = OLS_Robust(y(Break_Dates(3,1)+1:T),z(Break_Dates(3,1)+1:T,:),x(Break_Dates(3,1)+1:T,:), 0);
end


if size(z,2)==1
    Table_t_stat_dates=abs([t_stat1, 1, Break_Dates(1,1);
        t_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        t_stat3, Break_Dates(2,1)+1, Break_Dates(3,1);
        t_stat4, Break_Dates(3,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
else
    Table_t_stat_dates=abs([F_stat1, 1, Break_Dates(1,1);
        F_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        F_stat3, Break_Dates(2,1)+1, Break_Dates(3,1);
        F_stat4, Break_Dates(3,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
end


Missing=Tpi-A(ind_smaller,1)-1;
Start_date=Table_t_stat_dates(ind_larger,2);
End_date=Table_t_stat_dates(ind_larger,3);
Length_reg=End_date-Start_date;
for i=1:Length_reg-Missing
    if isempty(x)==1
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),[], 0);
    else
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),x(Start_date+i:Start_date + i +Missing,:), 0);
    end
end
    [value, ind]=min(vec_RSS);
    pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
    pi_Dates=sortrows(pi_Dates, 1,'ascend');
    pi_Dates_m3=pi_Dates;
    S_m3=pi_Dates_m3;
vec_RSS=zeros(T,1);


%%=========================================================================
%                               m=4 
%==========================================================================
Break_Dates=datevec(1:4,4);
if isempty(x)==1
[estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),[], 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),[], 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:Break_Dates(3,1)),z(Break_Dates(2,1)+1:Break_Dates(3,1),:),[], 0);
[estimator4, t_stat4, Fstat4] = OLS_Robust(y(Break_Dates(3,1)+1:Break_Dates(4,1)),z(Break_Dates(3,1)+1:Break_Dates(4,1),:),[], 0);
[estimator5, t_stat5, Fstat5] = OLS_Robust(y(Break_Dates(4,1)+1:T),z(Break_Dates(4,1)+1:T,:),[], 0);
else
    [estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),x(1:Break_Dates(1,1),:), 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),x(Break_Dates(1,1)+1:Break_Dates(2,1),:), 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:Break_Dates(3,1)),z(Break_Dates(2,1)+1:Break_Dates(3,1),:),x(Break_Dates(2,1)+1:Break_Dates(3,1),:), 0);
[estimator4, t_stat4, Fstat4] = OLS_Robust(y(Break_Dates(3,1)+1:T),z(Break_Dates(4,1)+1:T,:),x(Break_Dates(4,1)+1:T,:), 0);
[estimator5, t_stat5, Fstat5] = OLS_Robust(y(Break_Dates(4,1)+1:T),z(Break_Dates(4,1)+1:T,:),x(Break_Dates(4,1)+1:T,:), 0);
end


if size(z,2)==1
    Table_t_stat_dates=abs([t_stat1, 1, Break_Dates(1,1);
        t_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        t_stat3, Break_Dates(2,1)+1, Break_Dates(3,1);
        t_stat4, Break_Dates(3,1)+1, Break_Dates(4,1); ...
        t_stat5, Break_Dates(4,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
else
    Table_t_stat_dates=abs([F_stat1, 1, Break_Dates(1,1);
        F_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        F_stat3, Break_Dates(2,1)+1, Break_Dates(3,1);
        F_stat4, Break_Dates(3,1)+1, Break_Dates(4,1); ...
        F_stat5, Break_Dates(4,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
end


Missing=Tpi-A(ind_smaller,1)-1;
Start_date=Table_t_stat_dates(ind_larger,2);
End_date=Table_t_stat_dates(ind_larger,3);
Length_reg=End_date-Start_date;
for i=1:Length_reg-Missing
    if isempty(x)==1
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),[], 0);
    else
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),x(Start_date+i:Start_date + i +Missing,:), 0);
    end
end
    [value, ind]=min(vec_RSS);
    pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
    pi_Dates=sortrows(pi_Dates, 1,'ascend');
    pi_Dates_m4=pi_Dates;
    S_m4=pi_Dates_m4;
vec_RSS=zeros(T,1);


%%=========================================================================
%                               m=5 
%==========================================================================
Break_Dates=datevec(1:5,5);
if isempty(x)==1
[estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),[], 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),[], 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:Break_Dates(3,1)),z(Break_Dates(2,1)+1:Break_Dates(3,1),:),[], 0);
[estimator4, t_stat4, Fstat4] = OLS_Robust(y(Break_Dates(3,1)+1:Break_Dates(4,1)),z(Break_Dates(3,1)+1:Break_Dates(4,1),:),[], 0);
[estimator5, t_stat5, Fstat5] = OLS_Robust(y(Break_Dates(4,1)+1:Break_Dates(5,1)),z(Break_Dates(4,1)+1:Break_Dates(5,1),:),[], 0);
[estimator6, t_stat6, Fstat6] = OLS_Robust(y(Break_Dates(5,1)+1:T),z(Break_Dates(5,1)+1:T,:),[], 0);
else
    [estimator1, t_stat1, Fstat1] = OLS_Robust(y(1:Break_Dates(1,1)),z(1:Break_Dates(1,1),:),x(1:Break_Dates(1,1),:), 0);
[estimator2, t_stat2, Fstat2] = OLS_Robust(y(Break_Dates(1,1)+1:Break_Dates(2,1)),z(Break_Dates(1,1)+1:Break_Dates(2,1),:),x(Break_Dates(1,1)+1:Break_Dates(2,1),:), 0);
[estimator3, t_stat3, Fstat3] = OLS_Robust(y(Break_Dates(2,1)+1:Break_Dates(3,1)),z(Break_Dates(2,1)+1:Break_Dates(3,1),:),x(Break_Dates(2,1)+1:Break_Dates(3,1),:), 0);
[estimator4, t_stat4, Fstat4] = OLS_Robust(y(Break_Dates(3,1)+1:T),z(Break_Dates(4,1)+1:T,:),x(Break_Dates(4,1)+1:T,:), 0);
[estimator5, t_stat5, Fstat5] = OLS_Robust(y(Break_Dates(4,1)+1:Break_Dates(5,1)),z(Break_Dates(4,1)+1:Break_Dates(5,1),:),x(Break_Dates(4,1)+1:Break_Dates(5,1),:), 0);
[estimator6, t_stat6, Fstat6] = OLS_Robust(y(Break_Dates(5,1)+1:T),z(Break_Dates(5,1)+1:T,:),x(Break_Dates(5,1)+1:T,:), 0);

end


if size(z,2)==1
    Table_t_stat_dates=abs([t_stat1, 1, Break_Dates(1,1);
        t_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        t_stat3, Break_Dates(2,1)+1, Break_Dates(3,1);
        t_stat4, Break_Dates(3,1)+1, Break_Dates(4,1); ...
        t_stat5, Break_Dates(4,1)+1, Break_Dates(5,1);
        t_stat6, Break_Dates(5,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
else
    Table_t_stat_dates=abs([F_stat1, 1, Break_Dates(1,1);
        F_stat2, Break_Dates(1,1)+1, Break_Dates(2,1);
        F_stat3, Break_Dates(2,1)+1, Break_Dates(3,1);
        F_stat4, Break_Dates(3,1)+1, Break_Dates(4,1); ...
        F_stat5, Break_Dates(4,1)+1, Break_Dates(5,1);
        F_stat6, Break_Dates(5,1)+1, T]);
Table_t_stat_dates=sortrows(Table_t_stat_dates, 1,'descend');
A=[Table_t_stat_dates(:,3) - Table_t_stat_dates(:,2)];
A = cumsum( A , 1 );
Table_t_stat_dates=[Table_t_stat_dates A];
AA = A>Tpi;
AA = AA & cumsum(AA,1)<2;
[ind_larger,~] = find(AA);
ind_smaller=ind_larger-1;
end


Missing=Tpi-A(ind_smaller,1)-1;
Start_date=Table_t_stat_dates(ind_larger,2);
End_date=Table_t_stat_dates(ind_larger,3);
Length_reg=End_date-Start_date;
for i=1:Length_reg-Missing
    if isempty(x)==1
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),[], 0);
    else
       [~,~,~, vec_RSS(i,1)]=OLS_Robust(y(Start_date+i:Start_date + i +Missing),z(Start_date+i:Start_date + i +Missing,:),x(Start_date+i:Start_date + i +Missing,:), 0);
    end
end
    [value, ind]=min(vec_RSS);
    pi_Dates=[Table_t_stat_dates(1:ind_smaller,2:3); Start_date+ind, Start_date+ind+Missing-1];
    pi_Dates=sortrows(pi_Dates, 1,'ascend');
    pi_Dates_m5=pi_Dates;
    S_m5=pi_Dates_m5;
vec_RSS=zeros(T,1);