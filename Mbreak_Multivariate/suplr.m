% ***********************************************************************

function retp=suplr(maty,matx,m,n,br,rbeta,bigt,brbeta,brv)

% Procedure that constructs the supf test against m breaks.
% The test is the LR test under the optimal partitions 
% computed by dating_MLE.
% Control variables: 

% if brv==1, against breaks in variance
% if brbeta==1, against breaks in coefficients
% if brv==1 and brbeta==1, allowing changes in both.

pmatx=pzbar(matx,m,br,bigt);
res=maty-pmatx*rbeta; 
res=reshapeg(res,rows(res)/n,n);
seg=[0;br;(rows(maty)/n)];

if brv==1 && brbeta==1
        ulr=0;
        k=1;
        while k<=m+1
            i=seg(k)+1;
            j=seg(k+1);
            vvar=res(i:j,:)'*res(i:j,:)/(j-i+1);
            ulr=ulr-((j-i+1)*n/2)*(log(2*pi)+1)-(j-i+1)/2*log(det(vvar));
            k=k+1;
        end
end
if brv==1 && brbeta~=1
        %now, change in the var only, iterative method is used for the
         %optimization, with the unrestricted estimate as the prior input.
         %Then the beta is updated and then var. The procedure is iterated
         %until convergence. Non-linear optimization algorithm could be used
         %in the place of the iterative procedure
        
        ulr=0;
        k=1;
        while k<=m+1
            i=seg(k)+1;
            j=seg(k+1);
            vvar=res(i:j,:)'*res(i:j,:)/(j-i+1);
            ulr=ulr-((j-i+1)*n/2)*(log(2*pi)+1)-(j-i+1)/2*log(det(vvar));
            k=k+1;
        end
end
 if brbeta==1 &&  brv~=1
       ulr=0;
       vvar=res(1:bigt,:)'*res(1:bigt,:)/(bigt);
       ulr=-(bigt*n/2)*(log(2*pi)+1)-(bigt)/2*log(det(vvar));
 end
%now estimate the restricted model
    b=invpd(matx'*matx)*matx'*maty;
    res=maty-matx*b;
    umat=reshapeg(res,rows(res)/n,n); 
    vvar=umat'*umat/bigt;
    vstar=vvar+1;
    itr=1;
    while maxc(abs(vecr(vvar-vstar)))>1e-6 && itr<1000
	    vstar=vvar; 
        ibigv=kron(eye(rows(maty)/n),invpd(vstar)); 
	    b=invpd(matx'*ibigv*matx)*matx'*ibigv*maty; 
        res=maty-matx*b;
        umat=reshapeg(res,rows(res)/n,n); 
        vvar=umat'*umat/bigt;
	    itr=itr+1;
	    if itr==1000
            disp(' warning: the FGLS fail, the result does not converge');
        end
    end
    rlr=-(bigt*n/2)*(log(2*pi)+1)-(bigt)/2*log(det(vvar)); 
retp=2*(ulr-rlr);

