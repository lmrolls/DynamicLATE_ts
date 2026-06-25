% ********************************************************************@
% This procedure estimates the breaks and the coefficients 
% using the restrictions, based on the iterative method

function [dx,rbeta]=est(maty,matz,n,m,bigt,trm,R,brbeta,brv)
%  Input:
%         y: dependent variable
%         z: independent variable
%         q: total number of regression coefficients in a single regime
%         bigt: sample size
%   Output:
%         dx: the break dates
%         delta: the coefficients

    q=cols(matz);
    h=round(trm*bigt);
    tstar=bigt;  
	[global_,datevec,bigvec]=dating_mle(maty,matz,n,h,m,bigt);
                                     %find the optimal partiton of  the sample without restrictions 
                                     %on the coefficients, serves as the starting value of the 
                                     %iterative procedure@
    br=datevec(:,m);                   % @break dates@
    [rbeta,rvv]=r_estim(maty,matz,bigt,n,m,datevec(:,m),R,brbeta,brv);
	                                 %@now  find optimal partitions under the estimated coefficients@
	diff=0;
	while diff~=-1
             bigvec2=residuals(maty,matz,rbeta,q,n,m);
	        [global_,datevec,bigvec]=dating_m2(bigvec2,h,m,n,bigt);
	         if     datevec(:,m)==br
	                diff=-1;
			        dx=br;
             else
                     br=datevec(:,m);
                    [rbeta,rvv]=r_estim(maty,matz,bigt,n,m,br,R,brbeta,brv);
		             diff=diff+1;               
             end
    end
