% Find the optimal partition allowing m breaks, unrestricted model

function [global_,datevec,bigvec]=dating_mle(maty,matz,n,h,m,bigt)

% Input: 
% maty: vector dependent variables, with the (i-1)*n+j th element 
%       corresponding to the ith obseration from the jth equation
%       included in the regression
% matz: matrix of independent variables, arranged in the same manner
%       as maty. The number of columns corresponds to the number of
%       basic parameters included in the regression
% h:    mimunum number of observations in each segment
% m:    number of breaks.
% bigt: number of obervations in the sample
% Output: 
% global:  The maximized likelhood function under optimal partitons.
% datevec: the estimated break dates.
% bigvec:the big vector containing the MML for all possible segments


datevec=zeros(m,m); 

% up-tiangular matrix contains the estimated break dates for break 
% numbers from one to m

optdat=zeros(bigt,m); 

% row index corresponds to the ending dates, column index corresponds 
% to the number of breaks permitted before the ending date the cell 
% contains the optimal final break date

optmle=zeros(bigt,m);      

% same as above, the cell contains the MML corresponding to that break
% sturcture

dvec=zeros(bigt,1);                        

% the index is the date after which we inserting the break point.The 
% cell contains the corresponding MML

global_=zeros(m,1);                            

% Global MML when i breaks are permitted

bigvec=zeros(bigt*(bigt+1)/2,1);        

% the vector that contains the MML corresponding to all the possible segments. 
% The MML corresponding to the segment starting at J and lasting for span 
% corresponds to the index T(j-1)-(j-1)(j-2)/2+span


i=1;
while i <= bigt-h+1
    [vecmle]=mlef(i,maty,matz,h,n,bigt);
    bigvec((i-1)*bigt+i-(i-1)*i/2:i*bigt-(i-1)*i/2,1)=vecmle(i:bigt,1);
    i=i+1;
end


% Section that applies the dynamic programming algorithm to look for the
% optimal breaks
% The first case is when m = 1. Here the dynamic programming algorithm is not
% needed and one call to the parti(.) procedure is enough.


if m == 1
    [mlemax,datx]=parti(1,h,bigt-h,bigt,bigvec,bigt);
    datevec(1,1)=datx;
    global_(1,1)=mlemax; 
else

    %when m > 1, a dynamic programming algorithm is used.
    %The first step is to obtain the optimal one break partitions for all
    %possible ending dates from 2h to T-mh+1.
    %The optimal dates are stored in a vector optdat.
    %The associated MML are stored in a vector optssr.
    
    j1=2*h;                                                                              
    %First loop. Looking for the optimal one break partitions for break 
    %dates between h and T-h. j1 is the last date of the segments.
    
    while j1 <= bigt 
        [mlemax,datx]=parti(1,h,j1-h,j1,bigvec,bigt);           
        optmle(j1,1)=mlemax;            %@ again no typo@                                                     
        optdat(j1,1)=datx;                                                           
        j1=j1+1;
    end
    global_(1,1)=optmle(bigt,1);    %@ no typo@
    datevec(1,1)=optdat(bigt,1);

    
    %When this is done the algorithm looks for optimal 2,3,... breaks
    %partitions. The index used is ib.
    

    ib=2;
    while ib <= m
        if ib == m                 
            %if we have reached the highest number of breaks considered, only 
            %one segment is considered, that which ends at the last date of 
            %the sample.
            
            jlast=bigt;
            jb=ib*h;          %date of the break@
            while jb <=jlast-h
                dvec(jb,1) = optmle(jb,ib-1)+bigvec((jb+1)*bigt-jb*(jb+1)/2,1);
                jb=jb+1;
            end
            optmle(jlast,ib)=minc(dvec(ib*h:jlast-h,1));
            optdat(jlast,ib)=(ib*h-1)+minindc(dvec(ib*h:jlast-h,1));

        else
            %if we have not reached the highest number of breaks considered, 
            %we need to loop over the last date of the segment, between (ib+1)*h and T.
            
            jlast=(ib+1)*h;
            while jlast <= bigt
                jb=ib*h;             %date of the break@
                while jb <=jlast-h
                    dvec(jb,1) = optmle(jb,ib-1)+bigvec(jb*bigt-jb*(jb-1)/2+jlast-jb,1);
                    jb=jb+1;
                end
                optmle(jlast,ib)=minc(dvec(ib*h:jlast-h,1));
                optdat(jlast,ib)=(ib*h-1)+minindc(dvec(ib*h:jlast-h,1));
                jlast=jlast+1;
            end
        end

        %At each time we loop the results with ib breaks are retrieved
        %and printed@

        datevec(ib,ib)=optdat(bigt,ib);
        i=1;
        while i <= ib-1
            xx=ib-i;
            datevec(xx,ib)=optdat(datevec(xx+1,ib),xx);
            i=i+1;
        end
        global_(ib,1)=optmle(bigt,ib);
        ib=ib+1;
    end

end         %closing the if for the case m >1@

global_=-global_;

