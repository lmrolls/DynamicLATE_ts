% ********************************************************************

% Procedure that computes confidence intervals for the break dates based on 
% the "shrinking shifts asymptotic framework.  

function bound=interval2(maty,matx,m,n,bigt,br,beta,vv,hetq,vauto,brv,brbeta,prewhit)

% now constructing the quantities needed for the limit distribution

bound=zeros(m,4);
diagx=pzbar(matx,m,br,bigt);

br=[br;(rows(maty)/n)];
res=maty-diagx*beta;
res=reshapeg(res,rows(res)/n,n); 
% get standarized residuals
j=1;
while j<=m+1
    if j==1
        res(1:br(j),:)=res(1:br(j),:)*invpd(sqrm(vv(1:j*n,:)));
    else
        res(br(j-1)+1:br(j),:)=res(br(j-1)+1:br(j),:)*invpd(sqrm(vv((j-1)*n+1:j*n,:)));
    end
j=j+1;
end
% xy(seqa(1,1,rows(res)),res);
j=1;
while j<=m 
% construct confidence inteval for each break
vvar1=vv((j-1)*n+1:j*n,:);
vvar2=vv((j)*n+1:(j+1)*n,:); 
A1=sqrm(vvar1)*invpd(vvar2)*(vvar2-vvar1)*invpd(sqrm(vvar1));
A2=sqrm(vvar2)*invpd(vvar1)*(vvar2-vvar1)*invpd(sqrm(vvar2));
Q1=0;
Q2=0;
Pi1=0;
Pi2=0;
if j==1
    i=1;
else
    i=br(j-1)+1;
end
if hetq==0
    if (brv==1) && (vauto==1) % change in variance and/or in the regression coefficients, 
                                %    allow serial correlation
                                 %and change in the distribution of the regressors
        %construct estimate of pi1
        tempmat=zeros(br(j)-i+1,cols(matx));
        k=1;
        while k<=br(j)-i+1
            tempmat(k,:)=(matx((i-1+k-1)*n+1:(i-1+k)*n,:)'*invpd(vvar2)*sqrm(vvar1)*(res(i-1+k,:)'))';
            k=k+1;
        end
        pi1=correct(tempmat,prewhit);   %with scaling
        %@construct estimate of omiga@
        tempmat=zeros(br(j)-i+1,n^2);
        k=1;
        while k<=br(j)-i+1
            tempmat(k,:)=vec(res(i-1+k,:)'*res(i-1+k,:)-eye(n))';
            k=k+1;
        end
        omiga1=correct(tempmat,prewhit);       %with scaling
        while i<=br(j)
            Q1=Q1+matx((i-1)*n+1:i*n,:)'*invpd(vvar2)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        if j==1
            Q1=Q1/br(j);
        else
            Q1=Q1/(br(j)-br(j-1));
        end
        i=br(j)+1;
        %construct estimate of pi2
        tempmat=zeros(br(j+1)-i+1,cols(matx));
        k=1;
        while k<=br(j+1)-i+1
            tempmat(k,:)=(matx((i-1+k-1)*n+1:(i-1+k)*n,:)'*invpd(vvar1)*sqrm(vvar2)*(res(i-1+k,:)'))';
            k=k+1;
        end
        pi2=correct(tempmat,prewhit);   %with scaling
        tempmat=zeros(br(j+1)-i+1,n^2);
        k=1;
        while k<=br(j+1)-i+1
            tempmat(k,:)=vec(res(i-1+k,:)'*res(i-1+k,:)-eye(n))';
            k=k+1;
        end
        omiga2=correct(tempmat,prewhit);       %with scaling
        while i<=br(j+1)
            Q2=Q2+matx((i-1)*n+1:i*n,:)'*invpd(vvar1)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        Q2=Q2/(br(j+1)-br(j));
        diffb=beta(j*cols(matx)+1:(j+1)*cols(matx))-beta((j-1)*cols(matx)+1:(j)*cols(matx));
        Delt1=(1/2)*sumc(diag(A1*A1))+diffb'*Q1*diffb;
        Delt2=(1/2)*sumc(diag(A2*A2))+diffb'*Q2*diffb;
        Gam1=sqrt((1/4)*vec(A1)'*omiga1*vec(A1)+diffb'*pi1*diffb);
        Gam2=sqrt((1/4)*vec(A2)'*omiga2*vec(A2)+diffb'*pi2*diffb);
    end
    if (brv==1) && (vauto==0)  %change in the variance and/or in the regression coefficients
        omiga1=0;
        while i<=br(j)
            omiga1=omiga1+vec(res(i,:)'*res(i,:)-eye(n))*vec(res(i,:)'*res(i,:)-eye(n))';
            Q1=Q1+matx((i-1)*n+1:i*n,:)'*invpd(vvar2)*matx((i-1)*n+1:i*n,:);
            Pi1=Pi1+matx((i-1)*n+1:i*n,:)'*invpd(vvar2)*vvar1*invpd(vvar2)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        if j==1
            omiga1=omiga1/br(j);
            Q1=Q1/br(j);
            Pi1=Pi1/br(j);

% /*print omiga1 "mmiga";
% print Q1 "q1";
% print Pi1 "pi1";
% */

        else
            omiga1=omiga1/(br(j)-br(j-1));
            Q1=Q1/(br(j)-br(j-1));
            Pi1=Pi1/(br(j)-br(j-1));
        end
        i=br(j)+1;
        omiga2=0;
        while i<=br(j+1)
            omiga2=omiga2+vec(res(i,:)'*res(i,:)-eye(n))*vec(res(i,:)'*res(i,:)-eye(n))';
            Pi2=Pi2+matx((i-1)*n+1:(i)*n,:)'*invpd(vvar1)*vvar2*invpd(vvar1)*matx((i-1)*n+1:(i)*n,:);
            Q2=Q2+matx((i-1)*n+1:i*n,:)'*invpd(vvar1)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        omiga2=omiga2/(br(j+1)-br(j));
        Q2=Q2/(br(j+1)-br(j));
        Pi2=Pi2/(br(j+1)-br(j));
% /* print omiga2 "mmiga";
% print Q2 "q1";
% print Pi2 "pi1";
% */


        diffb=beta(j*cols(matx)+1:(j+1)*cols(matx))-beta((j-1)*cols(matx)+1:(j)*cols(matx));
        Delt1=(1/2)*sumc(diag(A1*A1))+diffb'*Q1*diffb;
        Delt2=(1/2)*sumc(diag(A2*A2))+diffb'*Q2*diffb;
        Gam1=sqrt((1/4)*vec(A1)'*omiga1*vec(A1)+diffb'*Pi1*diffb);
        Gam2=sqrt((1/4)*vec(A2)'*omiga2*vec(A2)+diffb'*Pi2*diffb);


    end
    if (brv==0) && (brbeta==1) && (vauto==1)  %change in the coefficients only

         %construct estimate of pi1
        tempmat=zeros(bigt,cols(matx));
        k=1;
        while k<=bigt
            tempmat(k,:)=(matx((k-1)*n+1:(k)*n,:)'*invpd(vvar2)*sqrm(vvar1)*(res(k,:)'))';
            k=k+1;
        end
        pi1=correct(tempmat,prewhit); %with scaling 
        k=1;
        while k<=bigt
            Q1=Q1+matx((k-1)*n+1:k*n,:)'*invpd(vvar2)*matx((k-1)*n+1:k*n,:);
            k=k+1;
        end
        Q1=Q1/bigt;
        pi2=pi1;
        Q2=Q1;
        diffb=beta(j*cols(matx)+1:(j+1)*cols(matx))-beta((j-1)*cols(matx)+1:(j)*cols(matx));
        Delt1=diffb'*Q1*diffb;
        Delt2=diffb'*Q2*diffb;
        Gam1=sqrt(diffb'*pi1*diffb);
        Gam2=sqrt(diffb'*pi2*diffb);
    end
    if (brv==0) && (brbeta==1) && (vauto==0)  %change in the coefficients only
        k=1;
        while k<=bigt
            Q1=Q1+matx((k-1)*n+1:k*n,:)'*invpd(vvar2)*matx((k-1)*n+1:k*n,:);
            Pi1=Pi1+matx((k-1)*n+1:k*n,:)'*invpd(vvar2)*vvar1*invpd(vvar2)*matx((k-1)*n+1:k*n,:);
            k=k+1;
        end
        Q1=Q1/bigt;
        Pi1=Pi1/bigt;
        Pi2=Pi1;
        Q2=Q1;
        diffb=beta(j*cols(matx)+1:(j+1)*cols(matx))-beta((j-1)*cols(matx)+1:(j)*cols(matx));
        Delt1=diffb'*Q1*diffb;
        Delt2=diffb'*Q2*diffb;
        Gam1=sqrt(diffb'*Pi1*diffb);
        Gam2=sqrt(diffb'*Pi2*diffb);
    end
elseif hetq==1     % allow the distribution of the regressors to be different, we only need to seperately treat vauto=1/0
                    % for a partial break model, the corresponding coefficients are zero
    if  vauto==1
        tempmat=zeros(br(j)-i+1,cols(matx)); 
        k=1;
        while k<=br(j)-i+1
            tempmat(k,:)=(matx((i-1+k-1)*n+1:(i-1+k)*n,:)'*invpd(vvar2)*sqrm(vvar1)*(res(i-1+k,:)'))';
            k=k+1;
        end
        pi1=correct(tempmat,prewhit);   %with scaling
        %construct estimate of omiga
        tempmat=zeros(br(j)-i+1,n^2);
        k=1;
        while k<=br(j)-i+1
            tempmat(k,:)=vec(res(i-1+k,:)'*res(i-1+k,:)-eye(n))';
            k=k+1;
        end
        omiga1=correct(tempmat,prewhit);       %with scaling     
        while i<=br(j)
            Q1=Q1+matx((i-1)*n+1:i*n,:)'*invpd(vvar2)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        if j==1
            Q1=Q1/br(j);
        else
            Q1=Q1/(br(j)-br(j-1));
        end
        i=br(j)+1;
        tempmat=zeros(br(j+1)-i+1,cols(matx));
        k=1;
        while k<=br(j+1)-i+1
            tempmat(k,:)=(matx((i-1+k-1)*n+1:(i-1+k)*n,:)'*invpd(vvar1)*sqrm(vvar2)*(res(i-1+k,:)'))';
            k=k+1;
        end
        pi2=correct(tempmat,prewhit); %with scaling
        tempmat=zeros(br(j+1)-i+1,n^2);
        k=1;
        while k<=br(j+1)-i+1
            tempmat(k,:)=vec(res(i-1+k,:)'*res(i-1+k,:)-eye(n))';
            k=k+1;
        end
        omiga2=correct(tempmat,prewhit);       %with scaling
        while i<=br(j+1)
            Q2=Q2+matx((i-1)*n+1:i*n,:)'*invpd(vvar1)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        Q2=Q2/(br(j+1)-br(j));
        diffb=beta(j*cols(matx)+1:(j+1)*cols(matx))-beta((j-1)*cols(matx)+1:(j)*cols(matx));
        Delt1=(1/2)*sumc(diag(A1*A1))+diffb'*Q1*diffb;
        Delt2=(1/2)*sumc(diag(A2*A2))+diffb'*Q2*diffb;
        Gam1=sqrt((1/4)*vec(A1)'*omiga1*vec(A1)+diffb'*pi1*diffb);
        Gam2=sqrt((1/4)*vec(A2)'*omiga2*vec(A2)+diffb'*pi2*diffb);
    end
    if vauto==0  %change in the variance and/or in the regression coefficients
        omiga1=0;
        to=0;
        while i<=br(j)
            omiga1=omiga1+vec(res(i,:)'*res(i,:)-eye(n))*vec(res(i,:)'*res(i,:)-eye(n))';
            Q1=Q1+matx((i-1)*n+1:i*n,:)'*invpd(vvar2)*matx((i-1)*n+1:i*n,:);
            Pi1=Pi1+matx((i-1)*n+1:i*n,:)'*invpd(vvar2)*vvar1*invpd(vvar2)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        if j==1
            omiga1=omiga1/br(j);
            Q1=Q1/br(j);
            Pi1=Pi1/br(j);
        else
            omiga1=omiga1/(br(j)-br(j-1));
            Q1=Q1/(br(j)-br(j-1));
            Pi1=Pi1/(br(j)-br(j-1));
        end
        i=br(j)+1;
        omiga2=0;
        while i<=br(j+1)
            omiga2=omiga2+vec(res(i,:)'*res(i,:)-eye(n))*vec(res(i,:)'*res(i,:)-eye(n))';
            Pi2=Pi2+matx((i-1)*n+1:(i)*n,:)'*invpd(vvar1)*vvar2*invpd(vvar1)*matx((i-1)*n+1:(i)*n,:);
            Q2=Q2+matx((i-1)*n+1:i*n,:)'*invpd(vvar1)*matx((i-1)*n+1:i*n,:);
            i=i+1;
        end
        omiga2=omiga2/(br(j+1)-br(j));
        Q2=Q2/(br(j+1)-br(j));
        Pi2=Pi2/(br(j+1)-br(j));
        diffb=beta(j*cols(matx)+1:(j+1)*cols(matx))-beta((j-1)*cols(matx)+1:(j)*cols(matx));
        Delt1=(1/2)*sumc(diag(A1*A1))+diffb'*Q1*diffb;
        Delt2=(1/2)*sumc(diag(A2*A2))+diffb'*Q2*diffb;
        Gam1=sqrt((1/4)*vec(A1)'*omiga1*vec(A1)+diffb'*Pi1*diffb);
        Gam2=sqrt((1/4)*vec(A2)'*omiga2*vec(A2)+diffb'*Pi2*diffb);
    end
end

eta=Delt2/Delt1;
phi1s=Gam1*Gam1/Delt1;
phi2s= Gam2*Gam2/Delt2;

cvf=cvg(eta,phi1s,phi2s);
a=(Delt1/Gam1)^2;
%print " a" a;
bound(j,1)=br(j,1)-cvf(4,1)/a;
bound(j,2)=br(j,1)-cvf(1,1)/a;
bound(j,3)=br(j,1)-cvf(3,1)/a;
bound(j,4)=br(j,1)-cvf(2,1)/a;
j=j+1;
end

bound(:,1)=floor(bound(:,1));
bound(:,2)=floor(bound(:,2))+1;
bound(:,3)=floor(bound(:,3));
bound(:,4)=floor(bound(:,4))+1;
