% ********************************************************************
%  confidence interval for the ordered breaks

function bound=interval3(matz,maty,n,start,last,dx1,dx2,delt1,delt2,res,segT,bigM,rep)

global hetq
step=1000;
obd=zeros(rep,3);
sigma=res'*res/rows(res);   %covariance matrix
if hetq==1
    ibigv=kron(eye(dx1-start+1),invpd(sigma)); 
    Qj=matz((start-1)*n+1:dx1*n,:)'*ibigv*matz((start-1)*n+1:dx1*n,:)/(dx1-start+1);
    ibigv=kron(eye(last-dx2),invpd(sigma)); 
    Qjp1=matz(dx2*n+1:last*n,:)'*ibigv*matz(dx2*n+1:last*n,:)/(last-dx2);
    P1j=delt1'*Qj*delt1;
    P2j=delt2'*Qj*delt2;
    P12j=delt1'*Qj*delt2;
    P1jp1=delt1'*Qjp1*delt1;
    P2jp1=delt2'*Qjp1*delt2;
    P12jp1=delt1'*Qjp1*delt2;
elseif hetq==0
    ibigv=kron(eye(last-start+1),invpd(sigma)); 
    Qj=matz((start-1)*n+1:last*n,:)'*ibigv*matz((start-1)*n+1:last*n,:)/(last-start);
    P1j=delt1'*Qj*delt1;
    P2j=delt2'*Qj*delt2;
    P12j=delt1'*Qj*delt2;
end
% print ";;" p1j p2j p12j p1jp1 p2jp1 p12jp1;
bigN=step*bigM;
interv=0.02; 
grids=1/interv;
gridN=grids*bigM;

simu=1;
while simu<=rep
    if find(simu==(10:10:rep))~=0
        disp(['rep: ' num2str(simu)]);
    end    
    bigW=rndn(bigN,4);
    bigW=cumsumc(bigW);
    bigW=bigW/sqrt(step);
    W=zeros(gridN,4);
    k=1;
    while k<=gridN
        W(k,:)=bigW(step*interv*(k-1)+1,:);
        k=k+1;
    end
    sv1=seqa(-bigM,1/grids,2*grids*bigM+1);
    sv2=sv1;
    V1=kron(sv1,ones(rows(sv1),1));
    V2=kron(ones(rows(sv2),1),sv2);
    sbv1=[rev(W(:,1));zeros(1,1);(W(:,4))];
    sbv2=[rev(W(:,2));zeros(1,1);(W(:,3))];  %two sided Brownian motion
    Bv1=kron(sbv1,ones(rows(sbv1),1));
    Bv2=kron(ones(rows(sbv2),1),sbv2);

    if hetq==1
        Z1=(V1<0).*(V2<=0).*(V1<=V2).*((Bv1+V1/2)+sqrt(P2j/P1j)*Bv2+V2/2*(P2j+2*P12j)/P1j);
        Z2=(V1<0).*(V2>0).*(Bv1+V1/2+sqrt(P2jp1/P1j)*Bv2-V2/2*(P2jp1/P1j));
        Z3=(V1>=0).*(V2>0).*(V1<=V2).*(sqrt(P1jp1/P1j)*Bv1+sqrt(P2jp1/P1j)*Bv2-(P2jp1/P1j)*V2/2-((P1jp1+2*P12jp1)/P1j)*V1/2);
        Z4=(V1>V2)*(-1);
    elseif hetq==0
        Z1=(V1<0).*(V2<=0).*(V1<=V2).*((Bv1+V1/2)+sqrt(P2j/P1j)*Bv2+V2/2*(P2j+2*P12j)/P1j);
        Z2=(V1<0).*(V2>0).*(Bv1+V1/2+sqrt(P2j/P1j)*Bv2-V2/2*(P2j/P1j));
        Z3=(V1>=0).*(V2>0).*(V1<=V2).*(Bv1+sqrt(P2j/P1j)*Bv2-(P2j/P1j)*v2/2-((P1j+2*P12j)/P1j)*V1/2);
        Z4=(V1>V2)*(-1);
    end
    Z=Z1+Z2+Z3+Z4;
    Z=[V1,V2,Z];
    Z=sortc(Z,3);
    obd(simu,:)=Z(rows(Z),:);
    simu=simu+1;

end
obd(:,1)=sortc(obd(:,1),1);
obd(:,2)=sortc(obd(:,2),1);%-(dx2-dx1)*p1j;
count=sumc(abs(abs(obd(:,1))-bigM)<1e-4)+sumc(abs(abs(obd(:,2))-bigM)<1e-4);
if count>=rep*0.05
    disp('wrong bound');
end
bound=zeros(4,2);
bound(1,:)=floor([dx1,dx2]-obd(floor(rep*0.025),1:2)/P1j)+ones(1,2);
bound(2,:)=floor([dx1,dx2]-obd(floor(rep*0.05),1:2)/P1j)+ones(1,2);
bound(3,:)=floor([dx1,dx2]-obd(floor(rep*0.95),1:2)/P1j);
bound(4,:)=floor([dx1,dx2]-obd(floor(rep*0.975),1:2)/P1j);

