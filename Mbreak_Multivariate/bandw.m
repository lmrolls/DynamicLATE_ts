% ********************************************************************@

function st=bandw(vhat)

% procedure that compute the automatic bandwidth based on AR(1)
% approximation for each vector of the matrix vhat. Each are given
% equal weigths of 1.@

nt=rows(vhat);
d=cols(vhat);
a2n=0;
a2d=0;
i=1;
while i <= d
    b=olsqr(vhat(2:nt,i),vhat(1:nt-1,i));
    sig=(vhat(2:nt,i)-b*vhat(1:nt-1,i))'*(vhat(2:nt,i)-b*vhat(1:nt-1,i));
    sig=sig/(nt-1);
    a2n=a2n+4*b*b*sig*sig/(1-b)^8;
    a2d=a2d+sig*sig/(1-b)^4;
    i=i+1;
end
a2=a2n/a2d;
st=1.3221*(a2*nt)^.2;
