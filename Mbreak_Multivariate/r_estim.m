% ***********************************************************************
% Procedrue that estimates the model allowing cross regime restrictions
% For the current version of the code, only linear restriction is allowed

function [nbeta,nvv]=r_estim(maty,matz,bigt,n,m,br,R,brbeta,brv)

% Input:
%         maty: dependent variable
%         matz: independent variable
%         n: number of equations
%         m: number of breaks
%         bigt: sample size
%         br: break dates
%         R: the restriction, taking the form beta=R*teta, where teta is the vector 
%         of basic parameters of the restricted model
%         brv: brv=1 if the variance-covariance matrix is allowed to change
%              brv=0 otherwise
%         brbeta: brbeta=1 if the coefficients are allowed to change
%                 brbeta=0 otherwise
%   Output:
%         nbeta: the estimate of coefficients imposing restrictions
%   Note: the model includes partial structural change and block-partial structural
%         change models
%         The estimation of based on iterative feasible GLS
 
seg=[0;br;(rows(maty)/n)];

% now estimate an unrestricted model

[beta,vv]=estim(maty,matz,n,m,br);
if brv==1 && brbeta==1
        pmatz=pzbar(matz,m,br,bigt);
        pmatz=pmatz*R;
        k=1;
        ibigv=eye(bigt*n);
        while k<=m+1
            i=seg(k)+1;
            j=seg(k+1);
            ibigv((i-1)*n+1:j*n,(i-1)*n+1:j*n)=kron(eye(j-i+1), invpd(vv((k-1)*n+1:k*n,:)));
            k=k+1;
        end
        b=invpd(pmatz'*ibigv*pmatz)*pmatz'*ibigv*maty;
        bstar=b+10;
        itr=1;
        while maxc(abs(bstar-b))>=1e-6 && itr<=1000
        %while maxc(abs(bstar-b))>=1 && itr<=1000
            bstar=b;
            k=1;
            %update the covariance matrix
            while k<=m+1
                i=seg(k)+1;
                j=seg(k+1);
                tempx=pmatz((i-1)*n+1:j*n,:);
	            tempy=maty((i-1)*n+1:j*n,:);
                res=tempy-tempx*bstar;
                umat=reshapeg(res,rows(res)/n,n); 
                vvar=umat'*umat/(j-i+1);
                vv((k-1)*n+1:k*n,:)=vvar;
                k=k+1;
            end
            k=1;
            % create the block diagonal matrix of the covariances
            while k<=m+1
                i=seg(k)+1;
                j=seg(k+1);
                ibigv((i-1)*n+1:j*n,(i-1)*n+1:j*n)=kron(eye(j-i+1),invpd(vv((k-1)*n+1:k*n,:)));
                k=k+1;
            end
            b=invpd(pmatz'*ibigv*pmatz)*pmatz'*ibigv*maty;
            itr=itr+1;
            if itr==1000
                disp(' The iteration has reached the upper limit');
            end
        end
    nbeta=R*b;
    nvv=vv;
end
        
if brv==1 && brbeta~=1
        %now, change in the covariance matrix only, iterative method is used for the
        % optimization, with the unrestricted estimate as the prior input.
        % The procedure is iterated until convergence. Non-linear optimization algorithm 
        % could be used in place of the iterative procedure
        
         ibigv=eye(bigt*n);
         k=1;
         while k<=m+1
            i=seg(k)+1;
            j=seg(k+1);
            ibigv((i-1)*n+1:j*n,(i-1)*n+1:j*n)=kron(eye(j-i+1),invpd(vv((k-1)*n+1:k*n,:)));
            k=k+1;
         end
         b=invpd(matz'*ibigv*matz)*matz'*ibigv*maty;
         bstar=b+10;
         itr=1;
         while maxc(abs(bstar-b))>=1e-6 && itr<=200
            bstar=b;
            %update the variance-covariance matrix
            k=1;
            while k<=m+1
                i=seg(k)+1;
                j=seg(k+1);
                tempx=matz((i-1)*n+1:j*n,:);
	            tempy=maty((i-1)*n+1:j*n,:);
                res=tempy-tempx*bstar;
                umat=reshapeg(res,rows(res)/n,n); 
                vvar=umat'*umat/(j-i+1);
                vv((k-1)*n+1:k*n,:)=vvar;
                k=k+1;
            end
            k=1;
            while k<=m+1
                i=seg(k)+1;
                j=seg(k+1);
                ibigv((i-1)*n+1:j*n,(i-1)*n+1:j*n)=kron(eye(j-i+1),invpd(vv((k-1)*n+1:k*n,:)));
                k=k+1;
            end
            b=invpd(matz'*ibigv*matz)*matz'*ibigv*maty;
            itr=itr+1;
            if itr==1000
                disp('The iteration has reached the upper limit');
            end
         end
    nbeta=b;
    nvv=vv;
    k=2;
    while k<=m+1
        nbeta=[nbeta;b];
        k=k+1;
    end
end

% the next allows changes in coefficients only
 if brbeta==1 &&  brv~=1
        k=1;
        vvar=0;
        while k<=m+1
            i=seg(k)+1;
            j=seg(k+1);
            tempx=matz((i-1)*n+1:j*n,:);
	        tempy=maty((i-1)*n+1:j*n,:);
            res=tempy-tempx*beta(k,:)';
            umat=reshapeg(res,rows(res)/n,n); 
            vvar=vvar+umat'*umat;
            k=k+1;
        end
            vvar=vvar/bigt;
            ibigv=kron(eye(bigt),invpd(vvar));
            pmatz=pzbar(matz,m,br,bigt);
            pmatz=pmatz*R;
            b=invpd(pmatz'*ibigv*pmatz)*pmatz'*ibigv*maty;
            %now iterate till convergence
            bstar=b+10;
        itr=1;
        while maxc(abs(vec(b-bstar)))>=1e-6 && itr<=1000
            bstar=b;
            k=1;
            res=maty-pmatz*bstar;
            umat=reshapeg(res,rows(res)/n,n); 
            vvar=umat'*umat/bigt;
            ibigv=kron(eye(bigt),invpd(vvar));
            b=invpd(pmatz'*ibigv*pmatz)*pmatz'*ibigv*maty;
            itr=itr+1;
            if itr==1000
                disp(' The iteration has reached the upper limit');
            end
        end
    nbeta=R*b;
    nvv=vvar;
    k=2;
    while k<=m+1
        nvv=[nvv;vvar];
        k=k+1;
    end
 end

