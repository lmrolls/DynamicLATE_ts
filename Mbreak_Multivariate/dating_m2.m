% @************************************************************************@
% /* This procedure is called when multiple breaks are present*/

function [global_,datevec,bigvec]=dating_m2(bigvec2,h,m,n,bigt)

% Input: 
% h:    mimunum number of observations in each segment
% m:    number of breaks.
% bigt: number of obervations in the sample
% Output: 
% global:  The maximized likelhood function under optimal partitons.
% datevec: the estimated break dates.

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
% @the vector that contains the MML corresponding to all the possible segments. 
% The MML corresponding to the segment starting at J and lasting for span 
% corresponds to the index T(j-1)-(j-1)(j-2)/2+span


if m == 1
    [mlemax,datx]=parti2(1,h,bigt-h,bigt,bigvec2,bigt,n);
    datevec(1,1)=datx;
    global_(1,1)=mlemax; 
else

    j1=2*h;                                                                              
    %First loop. Looking for the optimal one break partitions for break 
    %dates between h and T-h. j1 is the last date of the segments.
    
    while j1 <= bigt 
        [mlemax,datx]=parti2(1,h,j1-h,j1,bigvec2(1:2*bigt,:),bigt,n);           
        optmle(j1,1)=mlemax;         % again no typo                                                     
        optdat(j1,1)=datx;                                                           
        j1=j1+1;
    end
    global_(1,1)=optmle(bigt,1);   % no typo
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
            jb=ib*h;          %date of the break
            while jb <=jlast-h
                pmle=bigvec2(bigt*m+jb+1:bigt*(m+1),:)'*bigvec2(bigt*m+jb+1:bigt*(m+1),:)/(bigt-jb);
                dvec(jb,1) = optmle(jb,ib-1)+((bigt-jb+1)*n/2)*(log(2*pi)+1)+((bigt-jb+1)/2)*log(det(pmle));
                jb=jb+1;
            end
            optmle(jlast,ib)=minc(dvec(ib*h:jlast-h,1));
            optdat(jlast,ib)=(ib*h-1)+minindc(dvec(ib*h:jlast-h,1));

        else
            %if we have not reached the highest number of breaks considered, 
            %we need to loop over the last date of the segment, between (ib+1)*h and T.
            
            jlast=(ib+1)*h;
            while jlast <= bigt
                jb=ib*h;             %date of the break
                while jb <=jlast-h
                    pmle=bigvec2(bigt*ib+jb+1:bigt*ib+jlast,:)'*bigvec2(bigt*ib+jb+1:bigt*ib+jlast,:)/(jlast-jb);
                    dvec(jb,1) = optmle(jb,ib-1)+((jlast-jb+1)*n/2)*(log(2*pi)+1)+((jlast-jb+1)/2)*log(det(pmle));
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

