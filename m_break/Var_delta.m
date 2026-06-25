function vdel = Var_delta(y, z, i, q, bigt, b, prewhit, robust, x, p, withb, hetdat, hetvar)

% NAME            : Var_delta.m
% Description     : Compute Var-Cov matrix of delta
% Usage           : Var_delat(y, z, i, q, bigt, b, prewhit, robust, x, p, withb, hetdat, hetvar)
% Input           : y - dependent variable                                       
% Input           : z - matrix of regressors of dimension  Txq
% Input           : i - # of breaks under alternative
% Input           : q - # of regressors
% Input           : bigt - sample size
% Input           : b - Break dates obtained from Global method
% Input           : prewhit - pre-whitening                                    
% Input           : robust -  Robust method for heterogeneity and autocorrelation 
%                   in the residuals, 0 otherwise. The method used is Andrews(1991) 
%                   automatic bandwidth with AR(1) approximation and the quadratic 
%                   kernel. 
% Input           : x - regressors for partial change model
% Input           : p - # of regressors in x
% Input           : hetdat - hetdat=1; % option for the construction of the F tests. Set to 1 if want to
%                    allow different moment matrices of the regressors across segments. 
%                    If hetdat=0, the same moment matrices are assumed for each segment 
%                    and estimated from the ful sample. 
%                    It is recommended to set hetdat=1 if p>0.  regressors for partial change model
% Input           : hetvar - hetvar=1;  option for the construction of the F tests. Set to 1 if want to allow 
%                   for the variance of the residuals to be different across segments.
%                   If hetvar=0, the variance of the residuals is assumed constant 
%                   across segments and constructed from the ful sample. This option 
%                   is not available when robust=1.  
% Output          : SupFtest
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
% procedure that compute the covariance matrix of the estimates delta.

ev = ones(i+1,1);                                                               % # of regimes under H1
zbar = Partition_zbar(z,i,b);                                                   % Partition z

if p==0                                                                         % No x
    delv = olsqr(y,zbar);
    res = y - zbar*delv;
    reg = zbar;
else
    delv = olsqr(y,[zbar x]);
    res=y-[zbar x]*delv;
    
    if withb==0
        reg=zbar-x*inv(x'*x)*x'*zbar;
    else
        reg=[x zbar];
    end
end

if robust==0                                                                    %  robust -  Robust method for heterogeneity and autocorrelation 
                                                                                %                   in the residuals, 0 otherwise. The method used is Andrews(1991) 
                                                                                %                   automatic bandwidth with AR(1) approximation and the quadratic 
                                                                                %                   kernel.  
                                                                                
    
    if p==0
        if hetdat == 1 && hetvar == 0                                           % only regressors very across regimes
            sig=res'*res/bigt;
            vdel = sig*inv(reg'*reg);
        end
        if hetdat==1 && hetvar==1                                               % Also residuals' variance change across regimes 
            sig = Psigma(res, b, q, i, bigt);
            vdel=kron(sig,eye(q))*inv(reg'*reg);                                % So now use knronecker delta
        end
        if hetdat==0 && hetvar==0                                               % Just compute the average of 2nd moment of regressors
            lambda = Partitioned_lambda(b, i, bigt);
            sig=res'*res/bigt;
            vdel=sig*inv(kron(lambda,(z'*z)));
        end
        if hetdat==0 && hetvar==1
            lambda = Partitioned_lambda(b,i,bigt);
            sig = Psigma(res,b,q,i,bigt);
            vdel = kron(sig,eye(q))*inv(kron(lambda,(z'*z)));
        end
        
    else
        
        if hetdat==0
            disp('the case hetdat=0 is not allowed.');
            disp('vdel is returned as zeros.');
            vdel=zeros(q*(i+1),q*(i+1));
        end
        
        if hetdat==1 && hetvar==0
            sig=res'*res/bigt;
            vdel=sig*inv(reg'*reg);
        end
        
        if hetdat==1 && hetvar==1
            wbar=Partition_zbar(reg,i,b);
            ww=wbar'*wbar;
            sig=Psigma(res,b,q,i,bigt);
            gg=zeros((i+1)*q+p*withb, (i+1)*q+p*withb);
            ie=1;
            while ie<=i+1
                gg=gg+sig(ie,ie)*ww((ie-1)*((i+1)*q+p*withb)+1:ie*((i+1)*q+...
                    p*withb),(ie-1)*((i+1)*q+p*withb)+1:ie*((i+1)*q+p*withb));
                ie=ie+1;
            end
            vdel=inv(reg'*reg)*gg*inv(reg'*reg);
        end
    end
    
else
%==========================================================================
%               robust hac matrix is computed        
%==========================================================================  
    if hetdat==0
        disp('the case hetdat=0 is not allowed');
        disp('vdel is returned as zeros.');
        vdel=zeros(q*(i+1),q*(i+1));
    end
    
    if p==0
        if hetvar==1
            hac=zeros((i+1)*q,(i+1)*q);
            vdel=zeros((i+1)*q,(i+1)*q);
            hac(1:q,1:q)=b(1,1)*correct(z(1:b(1,1),:),res(1:b(1,1),1),prewhit);
            for j=2:i
                hac((j-1)*q+1:j*q,(j-1)*q+1:j*q)=(b(j,1)-b(j-1,1))*correct(z(b(j-1,1)+1:b(j,1),:),...
                    res(b(j-1,1)+1:b(j,1),1),prewhit);
            end
            
            hac(i*q+1:(i+1)*q,i*q+1:(i+1)*q)=(bigt-b(i,1))*correct(z(b(i,1)+1:bigt,:),...
                res(b(i,1)+1:bigt,1), prewhit);
            
            vdel=inv(reg'*reg)*hac*inv(reg'*reg);
        else
            hac=correct(z,res,prewhit);
            lambda=plambda(b,i,bigt);
            vdel=bigt*inv(reg'*reg)*kron(lambda,hac)*inv(reg'*reg);
        end
        
    else
        
         hac=correct(reg,res,prewhit);
         vdel=bigt*inv(reg'*reg)*hac*inv(reg'*reg);

        
    end
end

        
            
            