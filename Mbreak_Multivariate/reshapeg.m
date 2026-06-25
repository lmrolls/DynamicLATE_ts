function bmat=reshapeg(amat,a,b)

[r,c]=size(amat);

if r==1 && c==1
    
    bmat=amat*ones(a,b);

elseif r*c >=a*b 
    
    bmat=zeros(a,b);
    
    zvec=amat';
    zvec=zvec(:);
    kk=0;
    for j=1:a
        for k=1:b
           kk=kk+1;
            bmat(j,k)=zvec(kk,1);            
        end
    end
elseif r*c < a*b
    
    bmat=zeros(a,b);
    
    dd=ceil((a*b)/(r*c));
    
    if dd >= 2
        for jk=1:(dd-1)
            amat=[amat;amat];
        end
    end
    zvec=amat';
    zvec=zvec(:);
    kk=0;
    for j=1:a
        for k=1:b
           kk=kk+1;
            bmat(j,k)=zvec(kk,1);            
        end
    end

end    