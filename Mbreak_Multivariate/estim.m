function [beta,vv]=estim(maty,matz,n,m,br)

% procedure that estimates the coefficients under the 
% optimal partitions@

beta=zeros(m+1,cols(matz));
vv=zeros((m+1)*n,n);
br=[0;br;(rows(maty)/n)];
k=1;
while k<=m+1
    i=br(k)+1;
    j=br(k+1);
    segx=matz((i-1)*n+1:j*n,:);
    segy=maty((i-1)*n+1:j*n,:);
    b=invpd(segx'*segx)*segx'*segy;
    res=segy-segx*b;
    umat=reshapeg(res,rows(res)/n,n); 
    vvar=umat'*umat/(j-i+1);
    vstar=vvar+1;
    itr=1;
    while maxc(abs(vecr(vvar-vstar)))>1e-6 && itr<200
	    vstar=vvar; 
        ibigv=kron(eye(rows(segy)/n),invpd(vstar)); 
	    b=invpd(segx'*ibigv*segx)*segx'*ibigv*segy; 
        res=segy-segx*b;
        umat=reshapeg(res,rows(res)/n,n); 
        vvar=umat'*umat/(j-i+1);
	    itr=itr+1;
	    if itr==200
            disp(' warning: the FGLS fail, the result does not converge');
        end
    end
    beta(k,:)=b';
    vv((k-1)*n+1:k*n,:)=vvar;
    k=k+1;
end
