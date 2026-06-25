% ************************************************************************
% * the procedure is called only when number of breaks is one */
%  
% obtain an optimal one break partitions for a segment that starts at 
%     date start and ends at date last,under the given coefficients. 


function [optmle,dx]=parti2(start,b1,b2,last,bigvec2,bigt,n)

% Input: 
%         start: begining of the segment considered
%         b1: first possible break date
%         b2: last possible break date
%         last: end of segment considered
% Output: optmle: minus the value of the log-likelihood function
%         under the optimal partitons
%         dx: the break dates
% Note:   bigvec must be a matrix of size 2*bigt by n

dvec=zeros(bigt,1);      
j=b1;
while j<=b2
        m1=bigvec2(start:j,:)'*bigvec2(start:j,:)/(j-start+1);
        m2=bigvec2(1*bigt+j+1:1*bigt+last,:)'*bigvec2(1*bigt+j+1:1*bigt+last,:)/(last-j+1);
        dvec(j,1)=((last-start+1)*n/2)*(log(2*pi)+1)+((j-start+1)/2)*log(det(m1))+((last-j+1)/2)*log(det(m2)); 
        j=j+1;
end
optmle=minc(dvec(b1:b2,1));
dx=(b1-1)+minindc(dvec(b1:b2,1));
