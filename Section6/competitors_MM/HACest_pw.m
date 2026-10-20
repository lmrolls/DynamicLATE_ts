function V = HACest_pw(mU,mV,cLag)
%The function computes the HAC of mU,mV, using prewhitening and recolouring
%given the number of cLag.
%cLag ==-1 initial value for Newey-West 1994 algorithm for computing the lag 
%This function is better when we have two distict moment equations
%mU = T x p matrix of moments
%mV = T x q matrix of moments
%cLag=-1 HAC_lag = 4*[(T/100)^{2/9}]
%    =-2 HAC_lag = NW selection method

    
    ct=size(mU,1); 
    cp=size(mU,2);
    cq=size(mV,2);
    % Pre-whitening using one lag (see NW 1994)
    [mUhat, imA]=VARfilter(mU,cp);
    
    if cp==cq && all(mU(:)==mV(:))
        mVhat=mUhat;
        imB  =imA;
    else
        [mVhat, imB]=VARfilter(mV,cq);
    end
    
    function [mMhat, imA]=VARfilter(M,k)
        X      = M(1:ct-1,:);
        Y      = M(2:ct,:);
        A      =(Y'*X)/(X'*X);
        [U, S, V]= svd(A);
        S_tilde= diag(min(0.97,diag(S)));
        A      = U*S_tilde*V';
        imA    = eye(k)-A;
        mMhat  = Y-X*A';
    end
    
    V1 = HACest_tr(mUhat,mVhat,cLag);
    
    V  = (imA\V1)/imB'; % recolouring
    
    
end
