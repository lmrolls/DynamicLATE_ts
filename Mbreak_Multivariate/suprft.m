% ********************************************************************
% restricted supF test allowing m breaks 

function ftest=suprft(maty,matx,m,n,br,rbeta,R,bigt,brbeta,brv,prewhit) 

% local pmatx,res,xbstar,ubstar,ostar,sigma,umat,bigsig,tempmat,k,vmat,vbeta,vdelta;
% local rmat,rsub,j,ftest;

pmatx=pzbar(matx,m,br,bigt);
res=maty-pmatx*rbeta; 
umat=reshapeg(res,rows(res)/n,n); 
sigma=umat'*umat/bigt;
bigsig=invpd(sqrm(kron(eye(bigt),sigma)));
pmatx=pmatx*R;
xbstar=bigsig*pmatx;
ubstar=bigsig*res;
ubstar=reshapeg(ubstar,rows(ubstar)/n,n);
tempmat=zeros(bigt,cols(xbstar));
k=1;
while k<=bigt
    tempmat(k,:)=(xbstar((k-1)*n+1:(k)*n,:)'*(ubstar(k,:))')';
    k=k+1;
end
vmat=correct(tempmat,prewhit);
vdelta=invpd(xbstar'*xbstar/bigt)*(vmat/bigt)*invpd(xbstar'*xbstar/bigt);
vbeta=R*vdelta*R';
rsub=zeros(m,m+1);
j=1;
while j <= m
    rsub(j,j)=-1;
    rsub(j,j+1)=1;
    j=j+1;
end
rmat=kron(rsub, eye(cols(matx)));
ftest=rbeta'*rmat'*pinv(rmat*vbeta*rmat')*rmat*rbeta;
ftest=(bigt-(m+1)*cols(matx))*ftest/(bigt);
