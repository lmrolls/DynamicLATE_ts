% ********************************************************************
% sequential f test

function maxf=seqf(maty,matz,n,dt,nseg,bigt,h,brbeta,brv,prewhit)

global trm

% procedure that computes the supF(l+1|l) test l=nseg-1. The l breaks
% under the null are taken from the global minimization (in dt).

ssr=zeros(nseg,1);
lrtest=zeros(nseg,1);
dv=zeros(nseg+1,1);
dv(2:nseg,1)=dt;
dv(nseg+1,1)=bigt;
ds=zeros(nseg,1);

in=0;
is=1;
while is <= nseg
    length=dv(is+1,1)-dv(is,1);
    if length >= 2*h
        starti=dv(is,1)+1;
        endi=dv(is+1);
        segy=maty(n*(starti-1)+1:n*endi);
        segz=matz(n*(starti-1)+1:n*endi,:);
        [dxi,rbetai]=est(segy,segz,n,1,endi-starti+1,trm,eye(cols(matz)*2),1,0);
        lrtest(is)=supft(segy,segz,1,n,dxi,rbetai,endi-starti+1,1,0,prewhit);
    else
        in=in+1;
        lrtest(is,1)=0.0;
    end
    is=is+1;
end
if in == nseg
    disp('Given the location of the breaks from the global optimization');
    disp(['with ' num2str(nseg-1) ' breaks there was no more place to insert ']);
    disp('an additional breaks that satisfy the minimal length requirement.');
end
maxf=maxc(lrtest(1:nseg,1));
