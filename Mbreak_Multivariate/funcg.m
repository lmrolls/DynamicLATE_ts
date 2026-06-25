% *****************************************************************
% *obtain the quantiles of the limiting distribution

function g= funcg(x,bet,alph,b,deld,gam)

if x <= 0
    xb=bet*sqrt(abs(x));


    if abs(xb) <= 30 
        g=-sqrt(-x/(2*pi))*exp(x/8)...
        -(bet/alph)*exp(-alph*x)*cdfn(-bet*sqrt(abs(x)))...
        +((2*bet*bet/alph)-2-x/2)*cdfn(-sqrt(abs(x))/2);
    else
        aa=log(bet/alph)-alph*x-xb^2/2-log(sqrt(2*pi))-log(xb);
        g=-sqrt(-x/(2*pi))*exp(x/8)-exp(aa)*cdfn(-sqrt(abs(x))/2)...
        +((2*bet*bet/alph)-2-x/2)*cdfn(-sqrt(abs(x))/2);
    end

else
    xb=deld*sqrt(x);
    %print xb bet b deld gam*x -deld*sqrt(x) x;
    if abs(xb) <= 30
        g=1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)...
        +(b*deld/gam)*exp(gam*x)*cdfn(-deld*sqrt(x))...
        +(2-b*b*x/2-2*deld*deld/gam)*cdfn(-b*sqrt(x)/2);
    else
        aa=log((b*deld/gam))+gam*x-xb^2/2-log(sqrt(2*pi))-log(xb);
        g=1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)...
        +exp(aa)...
        +(2-b*b*x/2-2*deld*deld/gam)*cdfn(-b*sqrt(x)/2);

    end
end
