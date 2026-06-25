% ************************************************************************
function bigvec2=residuals(maty,matz,nbeta,q,n,m)

%  input:
%     y: dependent variable.
%     z: regressors
%     b: estimated coefficients.
%     q: total number of regression coefficients in each regime.
%     n: number of equations
%     m: number of breaks
%  output:
%     bigvec2: residuals for all the possible segments

bigt=rows(maty)/n;
bigvec2=zeros(bigt*(m+1),n);
i=1;
while i<=m+1
    bigvec2((i-1)*bigt+1:i*bigt,:)=reshapeg((maty-matz*nbeta((i-1)*q+1:i*q)),bigt,n);    
    i=i+1;
end
