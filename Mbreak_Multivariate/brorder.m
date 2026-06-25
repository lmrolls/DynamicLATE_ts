% ********************************************************************
%  Testing for ordered breaks, also estimating the break dates 

function [brdate,beta1,beta2,umat,lr]=brorder(maty,matz,n,q1,q2,h,T,trun)

TT=T*T;
vecmle=ones(TT,3+2*(q1+q2))*(-1e12);
vecmle2=ones(TT,3+2*(q1+q2))*(-1e12);
if cols(matz)~=q1+q2
    disp('wrong input of columns of regressors');
end
zz=[matz, matz];
z1=zz;
ii=h;
segx=zz;
while ii<=T-h
    ii=ii+1;
    jj=0;
    while (ii+jj)<=T-h
    z1(ii*n+1:rows(zz),cols(matz)+1:cols(matz)+q1)=zz(ii*n+1:rows(zz),cols(matz)+1:cols(matz)+q1);
    z1(1:ii*n,cols(matz)+1:cols(matz)+q1)=0.*zz(1:ii*n,cols(matz)+1:cols(matz)+q1);
    z1((ii+jj)*n+1:rows(zz),cols(matz)+q1+1:2*cols(matz))=zz((ii+jj)*n+1:rows(zz),cols(matz)+q1+1:2*cols(matz));
    z1(1:(ii+jj)*n,cols(matz)+q1+1:2*cols(matz))=0.*zz(1:(ii+jj)*n,cols(matz)+q1+1:2*cols(matz));
    segx=z1;
    segy=maty;
    if jj==0
        b=inv(segx'*segx)*segx'*segy;   
    else
        ibigv=kron(eye(rows(segy)/n),invpd(vvar)); 
        b=invpd(segx'*ibigv*segx)*segx'*ibigv*segy;
    end
    res=segy-segx*b;
    umat=reshapeg(res,rows(res)/n,n); 
    vvar=umat'*umat/rows(umat);
    vstar=vvar+1;
   itr=1;
    while maxc(abs(vecr(vvar-vstar)))>1e-8 && itr<1000
	    vstar=vvar; 
        ibigv=kron(eye(rows(segy)/n),invpd(vvar)); 
	    b=invpd(segx'*ibigv*segx)*segx'*ibigv*segy; 
        res=segy-segx*b;
        umat=reshapeg(res,rows(res)/n,n); 
        vvar=umat'*umat/rows(umat);
	    itr=itr+1; 
	    if itr==1000
	        disp(' warning: the FGLS fail, the result does not converge');
        end
    end 
    % now caculate the likelihood function for the given partition
    vecmle((ii-1)*T+jj,1)=ii;
    vecmle((ii-1)*T+jj,2)=jj+ii;
    vecmle((ii-1)*T+jj,3)=-T/2*log(det(vvar)); 
    vecmle((ii-1)*T+jj,3+1:3+2*(q1+q2))=b';
    jj=jj+1;
    end
end
vecmle2=vecmle;
vecmle2(:,3)=((vecmle2(:,2)-vecmle2(:,1))<=trun).*vecmle2(:,3)+((vecmle2(:,2)-vecmle2(:,1))>trun).*(-1e12);
vecmle=sortc(vecmle,3);
vecmle2=sortc(vecmle2,3);
H1like=vecmle2(rows(vecmle2),3);
% estimate a model with no change
    b=invpd(matz'*matz)*matz'*maty;  
    res=maty-matz*b; 
    umat=reshapeg(res,rows(res)/n,n); 
    vvar=umat'*umat/rows(umat);
    vstar=vvar+1;
   itr=1;
    while maxc(abs(vecr(vvar-vstar)))>1e-8 && itr<1000
	    vstar=vvar; 
        ibigv=kron(eye(rows(maty)/n),invpd(vvar)); 
	    b=invpd(matz'*ibigv*matz)*matz'*ibigv*maty; 
        res=maty-matz*b;
        umat=reshapeg(res,rows(res)/n,n); 
        vvar=umat'*umat/rows(umat);
	    itr=itr+1; 
	    if itr==1000
	        disp(' warning: the FGLS fail, the result does not converge');
        end
    end 
H0like=-T/2*log(det(vvar));
lr=2*(H1like-H0like);
brdate=vecmle(rows(vecmle),1:2);
beta1=vecmle(rows(vecmle),4:3+(q1+q2));
beta2=vecmle(rows(vecmle),4+(q1+q2):3+2*(q1+q2))+vecmle(rows(vecmle),4:3+(q1+q2));


beta1=beta1';
beta2=beta2';

