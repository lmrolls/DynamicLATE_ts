function []=mainop(y,z,n,q1,q2,trmm,bigt,hetq,trun,getcon,bigM,rep)

global T

h=floor(trmm*bigt);
disp('------------------------------------------------------------------------');
disp('The basic specifications for testing and estimation:');
disp('-------------------------------------------');
disp(['(1) Trimming= ' num2str(trmm)]);
disp(['(2) T= ' num2str(bigt)]);
disp('(3) The break in the first equation occurs first');
if hetq==1
    disp('(4) the distribution of the regressors is allowed to change');
else
    disp('(4) the distributions of the regressors is NOT allowed to change');
end
% transform the data
[maty,matz]=transf(y,z);
% testing and estimation
[brdate,beta1,beta2,res,lrtest]=brorder(maty,matz,n,q1,q2,h,bigt,trun);
delt1=[(beta2(1:q1)-beta1(1:q1));zeros(q2,1)];
delt2=[zeros(q1,1);(beta2(q1+1:q1+q2)-beta1(q1+1:q1+q2))];
if getcon==1
    bound=interval3(matz,maty,n,1,T,brdate(1),brdate(2),delt1,delt2,res,bigt,bigM,rep);
end
disp(' ');
disp('Output from the testing procedures:');
disp('------------------------------------------------------------------------');
clevel=[10,5,2.5,1];   %significance level 
cvs=cvlr(1,q1+q2,h/bigt);
disp(['The SupLR test for ordered breaks is ' num2str(lrtest)]);
disp('-----------------');
i=1;
while i<=4
    disp(['The critical values at the ' num2str(clevel(1,i)) '% level are']);
    disp(cvs(i));
    i=i+1;
end
disp('------------------------------------------------------------------------');
disp('Output from the estimation procedure');
disp('--------------------------------------------');
disp([' The estimated breaks are: ' num2str(brdate)]);
if getcon==1
    disp('--------------------------------------------');
    disp('Confidence intervals for the break dates are:');
    disp(['The 95% C.I. for the first break is: ' num2str([bound(4,1),bound(1,1)])]);
    disp(['The 95% C.I. for the second break is: ' num2str([bound(4,2),bound(1,2)])]);
    disp(['The 90% C.I. for the first break is: ' num2str([bound(3,1),bound(2,1)])]);
    disp(['The 90% C.I. for the second break is: ' num2str([bound(3,2),bound(2,2)])]);
end
disp('-------------------------------------------');
disp('The estimated coefficients are:');
disp(['For the first regime: ' num2str(beta1')]);
disp(['For the second regime: ' num2str(beta2')]);
disp('-------------------------------------------');
