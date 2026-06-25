function []=mainp(m,sumq,z,y,n,trm,bigt,brv,brbeta,vauto,hetq,prewhit)

global R 
% ***Input: m: number of breaks allowed
%           sumq: number of regression coefficients in each regime
%           z: regressors as a matrix, with the jth row containing 
%           observations at time j
%           y, dependent variable as a matrix, with the jth row containing 
%           observations at time j
%           n: the dimension of the dependent variable
%           trm: trimming 
%           bigT: sample size
%           brv: equal to 1 if covariance matrix is allowed to change
%           brbeta: equal to 1 if regression coefficients is allowed to change
%           vauto: equal to 1 if the error is allowed to be autocorrelated
%           hetq: equal to 1 if distribution of the regressors is allowed to change
%           prewhit: equal to 1 if use prewhittening

% ***Global: _S, matrix of constants selecting variables appearing in each 
%           regression.
%             R: matrix of constant imposing restrictions across regimes

disp('------------------------------------------------------------------------');
disp('The basic specifications for testing and estimation:');
disp('-------------------------------------------');
disp(['(1) M='  num2str(m)]); 

disp(['(2) Trimming=' num2str(trm)]); 

disp(['(3) T=' num2str(bigt)]); 

if brv==1
    disp('(4) The covariance matrix of errors is allowed to change');
    disp('    Normality is assumed when testing changes in the covariance matrix');
else
    disp('(4) The covariance matrix of errors is not allowed to change');
end

disp(['(5) The number of coefficients (beta) in each regime is: ' num2str(sumq)]); 

disp('-------------------------------------------');
disp('other specificiations:');
disp('-------------------------------------------');
disp('(1) the code reports results for an unrestricted model (apart');
disp('    from the  basic specifications and for a restricted model');
disp('    with restrictions specified by the user)');
if vauto==1
    disp('(2) the error is allowed to be autocorrelated');
else
    disp('(2) the error is serially uncorrelated');
end
if hetq==1
    disp('(3) the distribution of the regressors is allowed to change');
else
    disp('(3) the distributions of the regressors is NOT allowed to change');
end
if prewhit==1
    disp('(4) Pre-whittening with VAR(1) when constructing confidence intervals');
else
    disp('(4) No pre-whittening when constructing confidence intervals');
end

%disp '***********************************************************************';
disp('------------------------------------------------------------------------');
disp('Note: the critical values reported here assume q (the total number of'); 
disp(' parameters that are allowed to change) is not greater than 10.');   
disp('  For q>10, critical values can be simulated by applying Theorem 5'); 
disp('  in Qu and Perron (2007).') ;

h=round(trm*bigt);  %@trimming@
% transform the data
[maty,matz]=transf(y,z);
disp(' ');
disp('***********************************************************************');
disp('*     OUTPUT FROM THE UNRESTRICTED MODEL: TESTING AND ESTIMATION      *');
disp('***********************************************************************');

% the following calls the procedure to estimate and disp various
% tests for the number of breaks. It also returns the UDmax and WDmax
% tests. 

disp('Output from the testing procedures:');
disp('------------------------------------------------------------------------');
lrtest=zeros(m,1);
cvs=zeros(m,4);
clevel=[10,5,2.5,1];   %@significance level@

% /* The variance-covariance matrix of the errors is not allowed to change*/
if brv==0 && brbeta==1 
    disp('a) SupLR test for structural changes in the regression coefficients');
    disp(['(Number of changing parameters ' num2str(sumq) ')']);
    disp('-----------------');
    dxi=zeros(m,m);

    if vauto==0         %@no serial correlation@
        i=1;
        while i <= m
            [dxi(1:i,i),rbetai]=est(maty,matz,n,i,bigt,trm,eye(cols(matz)*(i+1)),1,0);
            lrtest(i)=suplr(maty,matz,i,n,dxi(1:i,i),rbetai,bigt,1,0);
            disp(['The SupLR test for 0 versus ' num2str(i) ' break ' num2str(lrtest(i,1))]); 
            i=i+1;
        end
        disp('-----------------');
        i=1;
        while i <= m
            cvs(i,:)=cvlr(i,sumq,trm);
            i=i+1;
        end
        i=1;
        while i<=4
            disp(['The critical values at the ' num2str(clevel(1,i)) '% level are (for 1 to ' num2str(m) ' breaks)']);
            disp(cvs(1:m,i)');
            i=i+1;
        end
        %The WD test
        disp('-------------------------------------------');
        disp(['b) The WDmax test for up to ' num2str(m) ' breaks']);
        %@disp '-------------------------------------------';@
        wd=zeros(1,4);
        i=1;
        while i<=m
            wd=(maxc([wd;(lrtest(i)*(cvs(1,:)./cvs(i,:)))]))';
            i=i+1;
        end
        cvs=cvwd(sumq,trm);
        i=1;
        while i<=4
            disp('-----------------');
            disp(['The WDmax test at the '  num2str(clevel(1,i)) ' %level is ' num2str(wd(i))]);
            disp(['(The critical value is ' num2str(cvs(i)) ' )']); 
            i=i+1;
        end
        disp('-------------------------------------------');
        disp('c) The output from the seq(l+1|l) test');
        i=1;
        while i<=m-1
            seq=seqlr(maty,matz,n,dxi(1:i,i),i+1,bigt,h,brbeta,brv);
            cvs=supseq(i,sumq,trm);
            disp('-----------------');
            disp(['The Seq(' num2str((i+1)) '|'  num2str(i) ') test is : ' num2str(seq)]); 
            disp(['(The critical values are ' num2str(cvs) ' )']);             
            i=i+1;
        end
    elseif vauto==1  %@Allowing serial correlation in the errors@
        i=1;
        while i <= m
            [dxi(1:i,i),rbetai]=est(maty,matz,n,i,bigt,trm,eye(cols(matz)*(i+1)),1,0);
            lrtest(i)=supft(maty,matz,i,n,dxi(1:i,i),rbetai,bigt,brbeta,brv,prewhit);
            disp(['The Supf test for 0 versus ' num2str(i) ' break ' num2str(lrtest(i,1))]); 
            i=i+1;
        end
        disp '-----------------';
        i=1;
        while i <= m
            cvs(i,:)=cvlr(i,sumq,trm);
            i=i+1;
        end
        i=1;
        while i<=4
            disp(['The critical values at the ' num2str(clevel(1,i)) '% level are (for 1 to ' num2str(m) ' breaks)']);
            disp(cvs(1:m,i)');
            i=i+1;
        end
        %@The WD test@
        disp('-------------------------------------------');
        disp(['b) The WDmax test for up to ' num2str(m) ' breaks']);
        %disp '-------------------------------------------';
        wd=zeros(1,4);
        i=1;
        while i<=m
            wd=(maxc([wd;(lrtest(i)*(cvs(1,:)./cvs(i,:)))]))';
            i=i+1;
        end
        cvs=cvwd(sumq,trm);
        i=1;
        while i<=4
            disp('-----------------');
            disp(['The WDmax test at the ' num2str(clevel(1,i)) ' %level is ' num2str(wd(i))]); 
            
            disp(['(The critical value is ' num2str(cvs(i)) ')']);
            i=i+1;
        end
        disp '-------------------------------------------';
        disp 'c) The output from the seq(l+1|l) test';
        i=1;
        while i<=m-1
            %@the sequential test@
            seq=seqf(maty,matz,n,dxi(1:i,i),i+1,bigt,h,brbeta,brv,prewhit);
            cvs=supseq(i,sumq,trm);
            disp('-----------------');
            disp(['The Seq(' num2str(i+1) '|' num2str(i)  ') test is : ' num2str(seq)]);
            disp(['(The critical values are '  num2str(cvs) ')']);
            i=i+1;
        end
    end
    disp '-------------------------------------------';
    disp '------------------------------------------------------------------------';
    disp(['Output from the estimation procedures allowing ' num2str(m) ' breaks']);
    [rbeta,rvv]=r_estim(maty,matz,bigt,n,m,dxi(:,m),eye(cols(matz)*(i+1)),brbeta,brv);
    bound2=interval2(maty,matz,m,n,bigt,dxi(:,m),rbeta,rvv,hetq,vauto,brv,brbeta,prewhit);
    disp('--------------------------------------------');
    disp([' The estimated breaks are: '  num2str(dxi(:,m)')]); 
    
    disp('--------------------------------------------');
    disp('Confidence intervals for the break dates are:');
    i=1;
    while i <= m
        disp(['The 95% C.I. for the ' num2str(i) 'th break is: ' num2str(bound2(i,1)) '  ' num2str(bound2(i,2))]); 
        
        disp(['The 90% C.I. for the ' num2str(i) 'th break is: ' num2str(bound2(i,3)) '  ' num2str(bound2(i,4))]); 
        
        i=i+1;
    end
    disp('--------------------------------------------');
    disp(' The estimated coefficients are:');
    jk=1;
    while jk<=m+1
        disp([' for the ' num2str(jk) 'th regime: '])  
        disp(num2str(rbeta((jk-1)*(rows(rbeta)/(m+1))+1:jk*(rows(rbeta)/(m+1)))'));
    jk=jk+1;
    end
    disp('---------------------------------------------');
    disp('The estimated covariance matrices of the errors are');
    jk=1;
    while jk<=m+1
        disp([' for the ' num2str(jk) 'th regime: ']) 
        disp(num2str(rvv((jk-1)*n+1:jk*n,:)));
        jk=jk+1;
    end
    disp '---------------------------------------------';
end

if brv==1 && brbeta==0
    lrtest=zeros(m,1);
    disp('-------------------------------------------');
    disp('a) SupLR test for pure structural change in covariance matrix (sigma)');
    disp(['Number of changing parameters ' num2str((n+1)*n/2)]);
    disp '-----------------';
    if vauto==1
        disp('NOTE! The code does not allow to test for structural change in the');
        disp('covariance matrix with serially correlated errors.');
        disp('The critical values are no longer valid');
        disp('-----------------');
    end
    i=1;
    dxi=zeros(m,m);
    while i <= m
    %@brbeta=0, brv=1@
        [dxi(1:i,i),rbetai]=est(maty,matz,n,i,bigt,trm,kron(ones((i+1),1),eye(cols(matz))),0,1);
        lrtest(i)=suplr(maty,matz,i,n,dxi(1:i,i),rbetai,bigt,0,1);
        disp(['The SupLR test for 0 versus '  num2str(i) ' break '  num2str(lrtest(i,1))]);
        i=i+1;
    end
    disp('-----------------');
    i=1;
    while i <= m
        cvs(i,:)=cvlr(i,(n+1)*n/2,trm);
        i=i+1;
    end
    i=1;
    while i<=4
        disp(['The critical values at the ' num2str(clevel(1,i)) ' %level are (for 0 against 1 to ' num2str(m) ' breaks)']);
        disp(num2str(cvs(1:m,i)'));
        i=i+1;
    end
    disp('-------------------------------------------');
    disp(['b) The WDmax test for up to ' num2str(m) ' breaks']);
    %@The WD test@
    wd=zeros(1,4);
    i=1;
    while i<=m
        wd=(maxc([wd;(lrtest(i)*(cvs(1,:)./cvs(i,:)))]))';
        i=i+1;
    end
    cvs=cvwd((n+1)*n/2,trm);
    i=1;
    while i<=4
        disp('-----------------');
        disp(['The WDmax test at the ' num2str(clevel(1,i)) ' %level is ' num2str(wd(i))]);
        disp(['(The critical value is ' num2str(cvs(i)) ')']);
        i=i+1;
    end
    disp('-------------------------------------------');
    disp('c) The output from the seq(l+1|l) test');
    i=1;
    while i<=m-1
        seq=seqlr(maty,matz,n,dxi(1:i,i),i+1,bigt,h,brbeta,brv);
        cvs=supseq(i,(n+1)*n/2,trm);
        disp('-----------------');
        disp(['The Seq(' num2str(i+1) '|' num2str(i) ') test is : '  num2str(seq)]);
        disp(['(The critical values are ' num2str(cvs) ')']);
        i=i+1;
    end
    disp('-------------------------------------------');
    disp('-------------------------------------------------------------------------');
    disp(['Output from the estimation procedures allowing ' num2str(m) ' breaks']);
    [rbeta,rvv]=r_estim(maty,matz,bigt,n,m,dxi(:,m),kron(ones((i+1),1),eye(cols(matz))),brbeta,brv);
    bound2=interval2(maty,matz,m,n,bigt,dxi(:,m),rbeta,rvv,hetq,vauto,brv,brbeta,prewhit);
    disp('--------------------------------------------');
    disp([' The estimated breaks are: ' num2str(dxi(:,m)')]); 
    disp('--------------------------------------------');
    disp('Confidence intervals for the break dates are:');
    i=1;
    while i <= m
        disp(['The 95% C.I. for the ' num2str(i) 'th break is: ' num2str(bound2(i,1)) '  ' num2str(bound2(i,2))]) 
        disp(['The 90% C.I. for the ' num2str(i) 'th break is: ' num2str(bound2(i,3)) '  ' num2str(bound2(i,4))]) 
        i=i+1;
    end
    disp('--------------------------- -----------------');
    disp(' The estimated coefficients are:');
    jk=1;
    while jk<=m+1
        disp([' for the ' num2str(jk) 'th regime: ']);  
        disp(num2str(rbeta((jk-1)*(rows(rbeta)/(m+1))+1:jk*(rows(rbeta)/(m+1)))'));
        jk=jk+1;
    end
    disp('---------------------------------------------');
    disp('The estimated covariance matrices of the errors are');
    jk=1;
    while jk<=m+1
        disp([' for the ' num2str(jk) 'th regime: ']); 
        disp(num2str(rvv((jk-1)*n+1:jk*n,:)));
        jk=jk+1;
    end
    disp('---------------------------------------------');
end

if brv==1 && brbeta==1
    lrtest=zeros(m,1);
    disp('-------------------------------------------');
    disp('a) SupLR test for pure structural change in regression coefficietns');
    disp(' (beta) and covariance matrix (sigma)');
    disp(['Number changing parameters' num2str(sumq+(n+1)*n/2)]);
    disp('-----------------');
    if vauto==1
        disp('NOTE: The code does not allow to test for structural change in the');
        disp('covariance matrix with serially correlated errors.');
        disp('The critical values are no longer valid');
        disp('-----------------');
    end
    dxi=zeros(m,m);
    i=1;
    while i <= m
    % brbeta=1, brv=1
    [dxi(1:i,i),rbetai]=est(maty,matz,n,i,bigt,trm,eye(cols(matz)*(i+1)),1,1);
    lrtest(i)=suplr(maty,matz,i,n,dxi(1:i,i),rbetai,bigt,1,1);
    disp(['The SupLR test for 0 versus ' num2str(i) ' break ' num2str(lrtest(i,1))]);
    i=i+1;
    end
    disp('-----------------');
    i=1;
    while i <= m
        cvs(i,:)=cvlr(i,sumq+(n+1)*n/2,trm);
        i=i+1;
    end
    i=1;
    while i<=4
        disp(['The critical values at the ' num2str(clevel(1,i)) '%level are (for 0 against 1 to ' num2str(m) ' breaks)']);
        disp(cvs(1:m,i)');
        i=i+1;
    end
    %@The WD test@
    disp('---------------------------------------------');
    disp(['b) The WDmax test for up to ' num2str(m) ' breaks']);
    wd=zeros(1,4);
    i=1;
    while i<=m
        wd=(maxc([wd;(lrtest(i)*(cvs(1,:)./cvs(i,:)))]))';
        i=i+1;
    end
    cvs=cvwd(sumq+(n+1)*n/2,trm);
    i=1;
    while i<=4
        disp('-----------------');
        disp(['The WDmax test at the ' num2str(clevel(1,i)) '% level is ' num2str(wd(i))]);
        disp(['(The critical value is ' num2str(cvs(i)) ')']);
        i=i+1;
    end
    disp('---------------------------------------------');
    disp('c) The output from the seq(l+1|l) test');
    i=1;
    while i<=m-1
        seq=seqlr(maty,matz,n,dxi(1:i,i),i+1,bigt,h,brbeta,brv);
        cvs=supseq(i,sumq+(n+1)*n/2,trm);
        disp('-----------------');
        disp(['The Seq( ' num2str(i+1) '|' num2str(i) ') test is : ' num2str(seq)]);
        disp(['(The critical values are ' num2str(cvs) ')']);
        i=i+1;
    end
    disp('-------------------------------------------');
    disp('-------------------------------------------------------------------------');
    disp(['Output from the estimation procedures allowing ' num2str(m) ' breaks']);
    [rbeta,rvv]=r_estim(maty,matz,bigt,n,m,dxi(:,m),eye(cols(matz)*(i+1)),brbeta,brv);
    bound2=interval2(maty,matz,m,n,bigt,dxi(:,m),rbeta,rvv,hetq,vauto,brv,brbeta,prewhit);
    disp('--------------------------------------------');
    disp([' The estimated breaks are: ' num2str(dxi(:,m)')]); 
    
    disp('--------------------------------------------');
    disp('Confidence intervals for the break dates are:');
    i=1;
    while i <= m
        disp(['The 95% C.I. for the ' num2str(i) 'th break is: ' num2str(bound2(i,1)) '  ' num2str(bound2(i,2))]);
        disp(['The 90% C.I. for the ' num2str(i) 'th break is: ' num2str(bound2(i,3)) '  ' num2str(bound2(i,4))]);
        i=i+1;
    end
    disp('--------------------------- -----------------');
    disp(' The estimated coefficients are:') ;
    jk=1;
    while jk<=m+1
        disp([' for the ' num2str(jk) ' th regime: '])  
        disp(rbeta((jk-1)*(rows(rbeta)/(m+1))+1:jk*(rows(rbeta)/(m+1)))');
        jk=jk+1;
    end
    disp('---------------------------------------------');
    disp('The estimated covariance matrices of the errors are');
    jk=1;
    while jk<=m+1
        disp([' for the ' num2str(jk) ' th regime: ']);
        disp(rvv((jk-1)*n+1:jk*n,:))
        jk=jk+1;
    end
    disp '---------------------------------------------';
end

disp(' ');
disp('***********************************************************************');
disp('*     OUTPUT FROM THE RESTRICTED MODEL: TESTING AND ESTIMATION      *');
disp('***********************************************************************');
disp('Output from the testing procedures:');
disp('------------------------------------------------------------------------');
[dx,rbeta]=est(maty,matz,n,m,bigt,trm,R,brbeta,brv);
[rbeta,rvv]=r_estim(maty,matz,bigt,n,m,dx,R,brbeta,brv);
if brv==1 && brbeta==0
    if vauto==0
        rlrtest=suplr(maty,matz,m,n,dx,rbeta,bigt,brbeta,brv);
        disp(['The total number coefficients allowed to change:' num2str(n*(n+1)/2)]);
        disp(['The value of the SUPF test: '  num2str(rlrtest)]);
        disp('-----------------');
        disp('The critical values are:');
        disp(num2str(cvlr(i,n*(n+1)/2,trm)));
    elseif vauto==1
        disp('NOTE: The code does not allow to test for structural change in the');
        disp('covariance matrix of the errors with serially correlated innovations.');
        disp('No result is reported');
        disp('-----------------');
    end
end
if brv==0 && brbeta==1
    cq=(n+1)*n/2*(brv==1)+rows(R)/(m+1)-(rows(R)-cols(R))/(m);
    disp(['The total number coefficients allowed to change: ' num2str(cq)]); 
    
    if vauto==0
        %no autocorrelation, the SUPLR test is used
        rlrtest=suplr(maty,matz,m,n,dx,rbeta,bigt,brbeta,brv);
        disp(['The value of the SUPLR test: ' num2str(rlrtest)]); 
        
        disp('-----------------');
    elseif vauto==1
        % use supf test with restrictions
        rftest=suprft(maty,matz,m,n,dx,rbeta,R,bigt,brbeta,brv,prewhit);
        disp(['The value of the SUPF test: '  num2str(rftest)]);
        disp('-----------------');
    end
    disp(['The critical values are: '  num2str(cvlr(i,cq,trm))]);
    
end
if brv==1 && brbeta==1
    if vauto==0
        cq=(n+1)*n/2*(brv==1)+rows(R)/(m+1)-(rows(R)-cols(R))/(m);
        disp(['The total number coefficients allowed to change: ' num2str(cq)]); 
        
        %@no autocorrelation, the SUPLR test is used@
        rlrtest=suplr(maty,matz,m,n,dx,rbeta,bigt,brbeta,brv);
        disp(['The value of the SUPLR test: ' num2str(rlrtest)]); 
        
        disp('-----------------');
        disp(['The critical values are: ' num2str(cvlr(i,cq,trm))]);
        
    elseif vauto==1
        disp('NOTE: The code does not allow to test for structural change in the');
        disp('covariance matrix of the errors with serially correlated innovations.');
        disp('Instead the test for change in regression coefficients is reported');
        disp('-----------------');
        cq=rows(R)/(m+1)-(rows(R)-cols(R))/(m);
        disp(['the total number of regression coefficients allowed to change: ' num2str(cq)]) 
        
        %@use supf test with restrictions@
        rftest=suprft(maty,matz,m,n,dx,rbeta,R,bigt,brbeta,brv,prewhit);
        disp(['The value of the SUPF test: ' num2str(rftest)]);
        disp('-----------------');
        disp(['The critical values are: ' num2str(cvlr(i,cq,trm))]);
        
    end
end
disp('------------------------------------------------------------------------');
disp(['Output from the estimation procedures allowing ' num2str(m) 'breaks']);
disp('             (The restricted model)             ');
disp('--------------------------------------------');
disp([' The estimated breaks are: ' num2str(dx')]); 

disp('--------------------------------------------');
bound2=interval2(maty,matz,m,n,bigt,dx,rbeta,rvv,hetq,vauto,brv,brbeta,prewhit);
disp('Confidence intervals for the break dates are:');
i=1;
while i <= m
    disp(['The 95% C.I. for the ' num2str(i) 'th break is: ' num2str([bound2(i,1),bound2(i,2)])]); 
    
    disp(['The 90% C.I. for the ' num2str(i) 'th break is: ' num2str([bound2(i,3), bound2(i,4)])]); 
    
    i=i+1;
end
disp('--------------------------- -----------------');
disp(' The estimated coefficients are:');
jk=1;
while jk<=m+1
    disp([' for the ' num2str(jk) 'th regime: '])  
    disp(rbeta((jk-1)*(rows(rbeta)/(m+1))+1:jk*(rows(rbeta)/(m+1)))');
jk=jk+1;
end
disp('---------------------------------------------');
disp('The estimated covariance matrices of the errors are');
jk=1;
while jk<=m+1
    disp([' for the ' num2str(jk) 'th regime: ']) 
    disp(rvv((jk-1)*n+1:jk*n,:));
    jk=jk+1;
end
disp('---------------------------------------------');
disp('***********************************************************************');
disp('*                              THE END!                                *');
disp('***********************************************************************');

