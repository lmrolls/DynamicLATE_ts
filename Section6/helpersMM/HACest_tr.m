function V = HACest_tr(mU,mV,cLag)
%The function computes the HAC of mU,mV, given the number of cLag.
%It is directly translated into MatLab from its Ox version.
%cLag ==-1 initial value for Newey-West 1994 algorithm for computing the lag 
%This function is better when we have two distict moment equations

    ct=size(mU,1); 
    cp=size(mU,2);
    cq=size(mV,2);
    
    %--- Save time in the next step
    if (cLag==-2)
        wu=ones(1,cp);
        wv=ones(1,cq);
    end
 
    if (size(mV,1)~=ct)
        return;
    end
    
    if (cLag<0)
        bauto=(cLag==-2);
        %--- O(T^(2/9))
        cLag=floor(4*(ct/100)^(2/9)); 
        if (bauto)
            %---When mV==mU, h = 2*mUwu. Since we use h to compute s0 and
            %s1, s1/s0 will cancel this extra2
            h = mU*wu'+mV*wv';
            s0=h'*h; 
            s1=0;
            for i=1:cLag
                sigi=h'*lag0(h,i);
                s0=s0+2*sigi;
                s1=s1+2*i*sigi;
            end
            %---bandwidth selection method
            cLag=floor(1.1447*(ct*(s1/s0)^2)^(1/3)); 
        end
    end
    
    if (cLag>=ct)
        cLag = ct-1; 
    end
    %---Bartlett Kernel - Newey and West (1987)
    hac=mU'*mV;
    for i=1:cLag
        kern=1-i/(cLag+1);                    
        hac=hac+kern*(lag0(mU,i)'*mV+mU'*lag0(mV,i));
    end
    
    V=hac/ct;
end
