% *************************************************************************@

function [ssrmin,dx]=parti(start,b1,b2,last,bigvec,bigt)

% procedure to obtain an optimal one break partitions for a segment that
% starts at date start and ends at date last. It return the optimal break
% and the associated minus likelhood value (MML).

% start: begining of the segment considered
% b1: first possible break date
% b2: last possible break date
% last: end of segment considered

dvec=zeros(bigt,1); 

% for a given ending date, the index is the inserting break point , 
% the cell contains the minus maximized likelhood function (MML)

ini=(start-1)*bigt-(start-2)*(start-1)/2+1;
%  the starting date of the data

j=b1;
while j<=b2
    jj=j-start;
    k=j*bigt-(j-1)*j/2+last-j;
    dvec(j,1)=bigvec(ini+jj,1)+bigvec(k,1);
    
    % the first term: the MML corresponding to data from date start to j
    % The second term: the ML corresponding to data from (j+1) and lasts 
    % for last -j
    
    j=j+1;
end
ssrmin=minc(dvec(b1:b2,1));
dx=(b1-1)+minindc(dvec(b1:b2,1));
