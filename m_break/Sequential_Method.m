function [nbreak,dv0] = Sequential_Method(m, signif,q,h,bigt,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)

% NAME            : Sequential_Method.m
% Description     : % This procedure applies the sequential procedure to obtain the number of
                    % breaks and the break dates. The current version only allows pure
                    % structural changes. This will be generalized.
                    % Here m is the maximum number of breaks allowed.      
% Usage           : Sequential_Method(m, signif,q,h,bigt,robust,prewhit,z,y,x,p,hetdat,hetvar,eps1)
% Input           : start - starting entry of the sample used.
% Input           : last - ending date of the last segment considered.
% Input           : y - dependent variable                                       
% Input           : z - matrix of regressors of dimension  Txq
% Input           : h - minimal length of a segment
% Output          : vec_Recursive_SSR - vector of sum of
%                   squared residuals (SSR) of length last-start+1 (stored for convenience in
%                   a vector of length T). 
%==========================================================================
% AUTHOR          : ALESSANDRO CASINI
% Email           : acasini@bu.edu  
% Date            : First version: 8/28/2015
% Last update     : 10/2/2015

%==========================================================================
% DISCLAIMER       
%==========================================================================  
% Please let me know if you find any bugs or encounter any problems while 
% using this code. All feedback is appreciated.
% @copyright, 2015

%==========================================================================
% REFERENCES
%==========================================================================
%  Bai, Jushan and Pierre Perron (1998): "Estimating and Testing Linear
%  Models with Multiple Structural Changes," Econometrica, vol 66, 
% and 
%  Bai, Jushan and Pierre Perron (2003): "Computation and Analysis of
%  Multiple Structural Change Models," Journal of Applied Econometrics, 18, 1-22.
% 
%
%==========================================================================
% DESCRIPTIONS: Computation Method
%==========================================================================
% this procedure computes recursive residuals from a data set that starts
% at date "start" and ends at date "last". It returns a vector of sum of
% squared residuals (SSR) of length last-start+1 (stored for convenience in
% a vector of length T).

% start: starting entry of the sample used.
% last: ending date of the last segment considered.
% y: dependent variable
% z: matrix of regressors of dimension Txq
% h: minimal length of a segment

                        
dv = zeros(m+2,1);
dv2 = zeros(m+2,1);
ftestv = zeros(m+1,1);

cv = getcv2(signif,eps1);                                                       % h = eps1*T
dv(1,1) = 0;

if p==0                                                                         % pure SC
    vssrev = Recursive_SSR(1,rot90(rot90(y)),rot90(rot90(z)), h, bigt);         % rot90 rotates a matrix 90 degrees
                                                                                % bigt=T;  % set effective sample size
                                                                                % This computes recursive residuals from data that starts at date 1 and end at date bigt=T
    vssr = Recursive_SSR(1,y,z,h,bigt);                                         % This computes recursive residuals from data that starts at date 1 and end at date bigt=T
    
    [ssrmin, datx]= Partition_1_Seq(h,bigt-h,bigt,vssr,vssrev);                 % procedure to obtain an optimal one break partitions for a segment that
                                                                                % starts at date start and ends at date last. It returns the optimal break
                                                                                % and the associated SSR.
                                                                                % Procedure used with the sequential method whenthe T*(T+1)/2 vector of SSR
                                                                                % is not computed.                     
else
    [ssrmin,datx]= onebp(y,z,x,h,1,bigt);
end

dv(2,1)=datx;

ftest = SupFTest_i_Breaks(y,z,1,q,bigt,dv(2,1),prewhit,robust,x,p,hetdat,hetvar);

if ftest < cv(q,1)
    nbreak=0;
    dv(2,1)=0;
    dv0 = dv(2:nbreak+1,1);
    nseg=1;
else
    disp(['The first break found is at: ' num2str(datx)]);
    nbreak=1;
    nseg=2;
    dv(nseg+1,1)=bigt;
end

while nseg <= m
    ds = zeros(nseg+1,1);
    ftestv = zeros(nseg+1,1);
    
    is=1;
    
    while is <= nseg
        length = dv(is+1,1) - dv(is,1);
        
        if length >= 2*h                                                        % lenght must be at least 2*h 
            
            if p==0
                vssr = Recursive_SSR(1,y(dv(is,1)+1:dv(is+1,1),1),...
                    z(dv(is,1)+1:dv(is+1,1),:),h,length);                       % This computes recursive residuals from data that starts at date 1=is+1 and end at date lenght
                
                vssrev = Recursive_SSR(1,rot90(rot90(y(dv(is,1)+1:dv(is+1,1),1))),...
                    rot90(rot90(z(dv(is,1)+1:dv(is+1,1),:))),h,length);         % This computes recursive residuals from data that starts at ROT date 1=is+1 and end at date lenght. These are the residuals for the 2nd half.
                
                [dum,ds(is,1)]=Partition_1_Seq(h,length-h,length,vssr,vssrev);  % procedure to obtain an optimal one break partitions for a segment that
                                                                                % starts at date start and ends at date last. It returns the optimal break
                                                                                % and the associated SSR.
                                                                                % Procedure used with the sequential method whenthe T*(T+1)/2 vector of SSR
                                                                                % is not computed.      
                ftestv(is,1) = SupFTest_i_Breaks(y(dv(is,1)+1:dv(is+1,1),1),...
                    z(dv(is,1)+1:dv(is+1,1),:),1,q,length,ds(is,1),...          % Sup test for 1 break
                    prewhit,robust,0,p,hetdat,hetvar);
            else
                [dum,ds(is,1)]=onebp(y,z,x,h,dv(is,1)+1,dv(is+1,1));
                ds(is,1)=ds(is,1)-dv(is,1);
                ftestv(is,1)=pftest(y(dv(is,1)+1:dv(is+1,1),1),...
                    z(dv(is,1)+1:dv(is+1,1),:),1,q,length,ds(is,1),...
                    prewhit,robust,x(dv(is,1)+1:dv(is+1,1),:),p,hetdat,hetvar);
            end
            
        else
            ftestv(is,1)=0.0;
        end
        
        is = is+1;
    end
    
    maxf = max(ftestv(1:nseg,1));                                               % Take the maximum ftestv
    
    if maxf < cv(q,nseg)
        nbreak = nbreak;
        dv0 = dv(2:nbreak+1,1);
    else
        [dum,newseg] = max(ftestv(1:nseg,1));
        disp(['The next break found is at: ' num2str(ds(newseg,1)+dv(newseg,1))]);
        
        dv(nseg+2,1)= ds(newseg,1)+dv(newseg,1);
        nbreak=nbreak+1;
        dv2 = sort(dv(2:nseg+2,1));
        dv(1,1)=0;
        dv(2:nseg+2,1)=dv2;
    end
    
    nseg=nseg+1;
end

disp(['The sequential procedure has reached the upper limit'])

dv0=dv(2:nbreak+1,1);
