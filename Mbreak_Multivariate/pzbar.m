% *****************************************************************************
%  procedure to construct the diagonal partition of matz with m break at dates b 
function matzb=pzbar(matz,m,bb,bigt)

nt=rows(matz);
q1=cols(matz);
n=nt/bigt;
if m==0
    disp(' m=0, no break is allowed, procedure terminates');
    matzb=0;
    return
else
    matzb=zeros(nt,(m+1)*q1);
    matzb(1:(bb(1,1)*n),1:q1)=matz(1:(bb(1,1)*n),:);
    i=2;
    while i <= m
        matzb(bb(i-1,1)*n+1:bb(i,1)*n,(i-1)*q1+1:i*q1)=matz(bb(i-1,1)*n+1:bb(i,1)*n,:);
        i=i+1;
    end
    matzb(bb(m,1)*n+1:nt,m*q1+1:(m+1)*q1)=matz(bb(m,1)*n+1:nt,:);
    
end

