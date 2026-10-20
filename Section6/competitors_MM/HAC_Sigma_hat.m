function Sigma_hat = HAC_Sigma_hat(Y_hat,Z,ky,kz,T)
% This function computes the HAC of vec(Z'*Yhat)
% Y_hat : Mx*Y matrix;
% Z     : instruments matrix;
% ky    : dimension of Y_hat
% kz    : dimension of Z
% T     : number of observations
% f_zv  : matrix of zeros(T,kz) auxiliary
% ons_kz: 1 x kz vector of ones auxiliary

pw   =  1; % 1 pre whitining for computing the full sample
cLag = -1; %-1 # of lags is proportional to sample size  [4*(T/100)^2] 
           %-2 # of lags is optimal according to Newey-West 
f_zv   = zeros(T,ky*kz);
ons_kz = ones(1,kz);          
           
h   = 1;
for j=1:ky
    f_zv(:,h:(j*kz)) = Z.*(Y_hat(:,j)*ons_kz);  
    h = j*kz+1;
end
f_zv      = f_zv -ones(T,1)*mean(f_zv);     % demeaned moments
Sigma_hat = HACest(f_zv,f_zv,cLag,pw);      %``prewhitening'' if pw==1
   
end
