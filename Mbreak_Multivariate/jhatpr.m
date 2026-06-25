function jhat=jhatpr(vmat)

% Procedure to compute the long-run covariance matrix of vmat.

nt=rows(vmat);
d=cols(vmat);
jhat=zeros(d,d);

% * calling the automatic bandwidth selection*/

st=bandw(vmat);

% lag 0 covariance
jhat=vmat'*vmat;

% forward sum
j=1;
while j <= nt-1
    jhat=jhat+kern(j/st)*vmat(j+1:nt,:)'*vmat(1:nt-j,:);
    j=j+1;
end

% bacward sum
j=1;
while j <= nt-1
    jhat=jhat+kern(j/st)*vmat(1:nt-j,:)'*vmat(j+1:nt,:);
    j=j+1;
end

% small sample correction@
jhat=jhat/(nt-d);

