function [glb datevec bigvec]= Dating_Breaks(y,z,h,m,q,bigt)

% NAME            : Dating_Breaks.m
% Description     : Compute the break points that globally minimizes the SSR.
%                   It returns optiomal dates and associated SSR for all
%                   numbers of breaks less than or equal to m. 
% Usage           : Dating_Breaks(y,z,h,m,q,bigt)
% Input           : start - beginning of the segment considered
% Input           : last - ending date of the last segment considered.
% Input           : b1 - first possible break date (usually set to h)                                    
% Input           : b2: last possible break date
% Output          : ssrmin - Minimum sum of squared residuals (SSR) 
% Output          : dx - Optimal break Date
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

                                                                                % This is the main procedure which calculates the break points that globally
                                                                                % minimizes the SSR. It returns optiomal dates and associated SSR for all
                                                                                % numbers of breaks less than or equal to m.
% =========================================================================
%             Prepare Space 
% =========================================================================
datevec = zeros(m,m);
optdat = zeros(bigt,m);
optssr = zeros(bigt,m);
dvec = zeros(bigt,1);
glb = zeros(m,1);
bigvec = zeros(bigt*(bigt+1)/2,1);


%% ========================================================================
%             Begin Algorithmn 
% =========================================================================
for i=1:bigt-h+1                                                                % Recall h is minimum size of each segment
    vecssr = Recursive_SSR(i,y,z,h,bigt);                                       % Call function reursive which for each i update the residuals
                                                                                % Computes recursive residuals from a data set that starts
                                                                                % at date "start" (i) and ends at date "last" (bigt). It returns a vector of sum of
                                                                                % squared residuals (SSR) of length last-start+1 (stored for convenience in
                                                                                % a vector of length T). 
    bigvec((i-1)*bigt+i-(i-1)*i/2:i*bigt-(i-1)*i/2,1) = vecssr(i:bigt,1);       % Store Recursive residuals in big vector
end

% =========================================================================
%            1 Break
% =========================================================================
if m==1                                                                         % If maximum number of breaks allowed is 1
    [ssrmin,datx] = Partition_1_Break(1, h, bigt-h, bigt, bigvec, bigt);        % Usage           : Partition_1 (start,b1,b2,last,bigvec,bigt)
                                                                                % Input           : start - beginning of the segment considered
                                                                                % Input           : last - ending date of the last segment considered.
                                                                                % Input           : b1 - first possible break date (usually set to h)                                    
                                                                                % Input           : b2: last possible break date
                                                                                % Input           : bigvec - vector storing recursive residuals (output
                                                                                %                   from Recursive_SSR.m
                                                                                % Input           : bigt - Sample size
                                                                                % Output          : ssrmin - Minimum sum of squared residuals (SSR) 
                                                                                % Output          : dx - Optimal break
                                                                                
    datevec(1,1) = datx;                                                        % Optimal break  
    glb(1,1) = ssrmin;                                                          % Minimum sum of squared residuals (SSR)              
    
else
% =========================================================================
%            >1 Break
% =========================================================================    
    for j1 = 2 * h:bigt                                                         % minimu break date is h. This computes optimal break date for each possible partitions from h to ...
                                                                                % But then we need to add residuals from pre-break.
        [ssrmin, datx] = Partition_1_Break(1, h, j1-h, j1, bigvec, bigt);       % Usage           : Partition_1 (start, b1, b2, last, bigvec, bigt)
                                                                                % Input           : 1 - beginning of the segment considered
                                                                                % Input           : last (bigt) - ending date of the last segment considered.
                                                                                % Input           : b1 (h) - first possible break date (usually set to h)                                    
                                                                                % Input           : b2 (j1-h) - last possible break date
                                                                                % Input           : bigvec - vector storing recursive residuals (output
                                                                                %                   from Recursive_SSR.m
                                                                                % Input           : bigt - Sample size
                                                                                % Output          : ssrmin - Minimum sum of squared residuals (SSR) 
                                                                                % Output          : datx - Optimal break
        optssr(j1,1)= ssrmin;
        optdat(j1,1)= datx;
    end
    glb(1,1)= optssr(bigt,1);                                                    % ssrmin - Minimum sum of squared residuals (SSR) BEFORE date bj                                                                                
    datevec(1,1)= optdat(bigt,1);                                                % datx - Optimal break BEFORE date bj
                                                                                 % Both stored in bigt at the end of procedure
    
    for ib=2:m
        %------------------------------------------------------------------
        %           For Max # of Breaks
        %------------------------------------------------------------------        
        if ib==m
            jlast = bigt;
            for jb = ib*h:jlast-h                                                 % jb starts at h*2 at least        
                dvec(jb,1) = optssr(jb,ib-1) + bigvec((jb+1)*bigt-jb*(jb+1)/2,1); % Add residuals from pre-break  
            end
            optssr(jlast,ib)=min(dvec(ib*h:jlast-h,1))';
            [dummy, minindcdvec]=min(dvec(ib*h:jlast-h,1));
            optdat(jlast,ib)=(ib*h-1)+minindcdvec';
            
        else
        %------------------------------------------------------------------
        %           For <Max # of Breaks
        %------------------------------------------------------------------            
            for jlast = (ib+1)*h:bigt                                             % Last Potential break dates to try after jb
                for jb = ib*h: jlast-h                                            % Potential break dates to try after jb
                    dvec(jb,1) = optssr(jb,ib-1) + bigvec(jb*bigt-jb*(jb-1)/2+jlast-jb,1); % SSR: 1st term is SSR pre-jb and 2nd term is after jb
                end
                optssr(jlast,ib) = min(dvec(ib*h:jlast-h,1));                   % For each possible jlast
                [dummy, minindcdvec] = min(dvec(ib*h:jlast-h,1));
                optdat(jlast,ib) = (ib*h-1) + minindcdvec';                     % Date of last break
            end
        end
       
        datevec(ib,ib) = optdat(bigt,ib);                                       % This picks last date of break for each m
        
        for i = 1:ib-1
            xx = ib-i;                                                            % xx goes from 1 to ib-1 
            datevec(xx,ib) = optdat(datevec(xx+1,ib),xx);
        end
        glb(ib,1) = optssr(bigt,ib);
    end
end

    