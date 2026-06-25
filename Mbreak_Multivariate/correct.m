% ********************************************************************@
% 
% * procedure that obtains the robust standard errors*/
function hac=correct(res,prewhit)

% main procedures which activates the computation of robust standard
% errors.

d=cols(res);
nt=rows(res);
b=zeros(d,1);
bmat=zeros(d,d);
vstar=zeros(nt-1,d);
vmat=res;

% Procedure that applies prewhitenning to the matrix vmat by filtering
% with a VAR(1). if prewhit=0, it is skipped.


if prewhit == 1

    % /* estimating the VAR(1) */
    i=1;
    while  i <= d
        b=olsqr(vmat(2:nt,i),vmat(1:nt-1,:));
        bmat(:,i)=b;
        vstar(:,i)=vmat(2:nt,i)-vmat(1:nt-1,:)*b;
        i=i+1;
    end

    % /* call the kernel on the residuals */

    jh=jhatpr(vstar);

    % /* recolor */

    hac=inv(eye(d)-bmat)*jh*(inv(eye(d)-bmat))';

else
    hac=jhatpr(vmat);
end

