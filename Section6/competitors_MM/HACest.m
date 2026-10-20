function V = HACest(mU,mV,cLag, pw)
%The function computes the HAC of mU,mV, 
%pw==1: use prewhitening and recolouring 
%given the number of cLag.
%cLag ==-1 initial value for Newey-West 1994 algorithm for computing the lag 
%This function is better when we have two distict moment equations
%mU = T x p matrix of moments
%mV = T x q matrix of moments
%cLag=-1 HAC_lag = 4*[(T/100)^{2/9}]
%    =-2 HAC_lag = NW selection method

    if (pw==1)
            V= HACest_pw(mU,mV,cLag);
    else
            V= HACest_tr(mU,mV,cLag);
    end
    
    
end
