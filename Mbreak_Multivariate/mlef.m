% ************************************************************************

% /*Computing MLE based on a recursive formula*/
function vecmle=mlef(start,maty,matz,h,n,last)

% ***Function: This procedure computes MLE from a data set that
%              starts at date 'start' and ends at date 'last'. 
% ***return:   It returns a vector of MLE of length last-start+1
%              (stored for convenience in a vector of length T).
% ***method:   iterated FGLS is used to obtain the estimates, the recursive formula
%              in Tobing and McGilchrist (1992) is applied to reduce
%              the computational cost.
% 
% ***Inputs:
%             start: starting entry of the sample used.
%             last: ending date of the last segment considered
%             maty: vector of dependent variables, with the n*(j-1)+i th element 
%                   corresponding to the jth observation from ith equation. 
%             matz: matrix of all regressors included in the regression, arranged 
%                   in the same way as maty
%             h:    minimal length of a segment
%             n:    number of equations included in the system.
% 
% ***Note: that the first h-1 elements have no meaning, usually is 
%          initialized as zero
% ***Note: the maximum number of iterations allowed is 1000


vecmle=zeros(last,1);
% initialize the beta and sigma using the first h observations, 
% iterated until convergence

i=start;
j=start+h-1;
segx=matz((i-1)*n+1:j*n,:);
segy=maty((i-1)*n+1:j*n,:);
b=invpd(segx'*segx)*segx'*segy;
res=segy-segx*b;
umat=reshapeg(res,rows(res)/n,n); 
vvar=umat'*umat/h;
vstar=vvar+1;
itr=1;
while maxc(abs(vecr(vvar-vstar)))>1e-6 && itr<1000
    %while maxc(abs(vecr(vvar-vstar)))>1 && itr<1000
	vstar=vvar; 
    ibigv=kron(eye(rows(segy)/n),invpd(vstar)); 
	b=invpd(segx'*ibigv*segx)*segx'*ibigv*segy; 
    res=segy-segx*b;
    umat=reshapeg(res,rows(res)/n,n); 
    vvar=umat'*umat/h;
	itr=itr+1;
	if itr==1000
        disp(' warning: The number of iterations has reached the limit, FGLS does not converge');
        disp(' consider to increase the minimum length of a segment');
    end
end
bstar=b;
ibigv=kron(eye(h),invpd(vvar));
ihstar=invpd(segx'*ibigv*segx);
bstar=b;
ystar=segy;
zstar=segx;
vstar=vvar;

% start recursion

while j<=last
	gstar=ihstar*matz((j-1)*n+1:j*n,:)';
	icstar=invpd(vstar+matz((j-1)*n+1:j*n,:)*ihstar*matz((j-1)*n+1:j*n,:)'); 
	b=bstar+gstar*icstar*(maty((j-1)*n+1:j*n,:)-matz((j-1)*n+1:j*n,:)*bstar);
	tempz=[zstar;matz((j-1)*n+1:j*n,:)];
	tempy=[ystar;maty((j-1)*n+1:j*n,:)];
    res=tempy-tempz*b;
    umat=reshapeg(res,rows(res)/n,n); 
    vvar=umat'*umat/(j-i+1);
	itr=1;
    vstar=vvar+10;
    while maxc(abs(vecr(vvar-vstar)))>1e-6 && itr<1000
	    vstar=vvar;
        
        % note that gstar does not need to be updated
	    
        icstar=invpd(vstar+matz((j-1)*n+1:j*n,:)*ihstar*matz((j-1)*n+1:j*n,:)');
	    b=bstar+gstar*icstar*(maty((j-1)*n+1:j*n,:)-matz((j-1)*n+1:j*n,:)*bstar);
        res=tempy-tempz*b;
        umat=reshapeg(res,rows(res)/n,n); 
        vvar=umat'*umat/(j-i+1);
	    itr=itr+1;
	    if itr==1000
            disp(' warning: The number of iterations has reached the limit, FGLS does not converge');
            disp(' consider to increase the minimum length of a segment');
        end
    end
    vecmle(j,1)= ((j-i+1)*n/2)*(log(2*pi)+1)+(j-i+1)/2*log(det(vvar)); 
    % note that minus the likelihood function is recorded
    bstar=b;
    zstar=tempz;
    ystar=tempy;
    vstar=vvar;
    ihstar=ihstar-gstar*icstar*gstar';
    j=j+1;
end

