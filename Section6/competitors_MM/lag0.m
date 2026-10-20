function lagv=lag0(v,i)
%This function returns the i_th lag of vector/matrix v

%Input
%v   : a vector or matrix
%i   : number of lag

%Output
%lagv: the i_th lag of v

    m=size(v,1);      % no. of row of v
    n=size(v,2);      % no. of column of v
    lagv=zeros(m,n);  % initialize lagv, with the same size as v
    %Note: this function is only called to compute inner products, thus it 
    %is ok to use "0" as missing values. 
    
    %replace lagv with its i_th lag
    lagv(i+1:m,:)=v(1:m-i,:); 
end