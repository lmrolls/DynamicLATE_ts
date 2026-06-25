% ********************************************************************
% *supF test allowing m breaks

function ftest=supft(maty,matx,m,n,br,rbeta,bigt,brbeta,brv,prewhit)

pmatx=pzbar(matx,m,br,bigt);
res=maty-pmatx*rbeta; 
umat=reshapeg(res,rows(res)/n,n); 
sigma=umat'*umat/bigt;
bigsig=invpd(sqrm(kron(eye(bigt),sigma)));
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
vbeta=invpd(xbstar'*xbstar/bigt)*(vmat/bigt)*invpd(xbstar'*xbstar/bigt);
rsub=zeros(m,m+1);
j=1;
while j <= m
    rsub(j,j)=-1;
    rsub(j,j+1)=1;
    j=j+1;
end
rmat=kron(rsub,eye(cols(matx)));
ftest=rbeta'*rmat'*inv(rmat*vbeta*rmat')*rmat*rbeta;
ftest=(bigt-(m+1)*cols(matx))*ftest/(bigt);
