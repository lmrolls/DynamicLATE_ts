function [glob,datevec,bigvec]= Procedure_piL_Full_Sample(bigt,y,z,q,m,h,eps1,robust,prewhit,hetomega,hetq,doglobal,fixb,x,p,eps,maxi,betaini,printd,hetdat,hetvar)
% 
% disp('The option chosen are:')
% disp(['h = ' num2str(h)])
% disp(['eps1 = ' num2str(eps1)])
% disp(['hetdat = ' num2str(hetdat)])
% disp(['hetvar = ' num2str(hetvar)])
% disp(['hetomega = ' num2str(hetomega)])
% disp(['hetq = ' num2str(hetq)])
% disp(['robust = ' num2str(robust) ' (prewhit =' num2str(prewhit) ')'])
% disp(['hetvar = ' num2str(hetvar)])
% disp(['The maximum number of breaks is: ' num2str(m)])
% disp('***********************************************************')

siglev=[10;5;2.5;1];

%% ========================================================================
%               Global Method       
%==========================================================================
if doglobal==1                                                                  % Call the global procedure to obtain the break dates and the associated sum of squared residuals for all numbers of breaks between 1 and m.                                                                    
    
 % disp('Output from the global optimization')
 % disp('***********************************************************')
 
% -------------------------------------------------------------------------
%              Pure Structural Change Model      
% -------------------------------------------------------------------------
 if p==0                                                                       % p=0 calls pure change model
    [glob, datevec, bigvec] = Dating_Breaks(y, z, h, m, q, bigt);                   
                                                                               % glob = global minimizer SSR
                                                                               % datevec = global minimum break dates 
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
 else
% -------------------------------------------------------------------------
%              Partial Structural Change Model      
% -------------------------------------------------------------------------
     % disp('This is a partial structural change model and the following')
     % disp('specifications were used:')
     % disp(['The number of regressors with fixed coefficients, p: ' num2str(p)])
     % disp(['Whether (1) or not (0) the initial values of beta are fixed: ' num2str(fixb)])
     % disp(['If so, at the following values: ' num2str(betaini)])
     % disp(['The convergence criterion is: ' num2str(eps)])
     % disp(['Whether (1) or not (0) the iterations are printed: ' num2str(printd)])
     % disp('----------------------')
     % 
     [glob,datevec,bigvec] = nldat(y, z, x, h, m, p, q, bigt, fixb, eps, maxi, betaini, printd);
 end
  
 for i=1:m
     % disp(['The model with ' num2str(i) 'breaks has SSR : ' num2str(glob(i,1))])
     % disp(['The dates of the breaks are: ' num2str(datevec(1:i,i)')])
 end
end
