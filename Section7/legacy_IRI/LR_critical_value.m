function cv_LR=LR_critical_value(M2_ST, q, alpha)
n_reps=1000;
for i_reps=1:n_reps
    M1_ST=chi2rnd(q);
    M12_ST=sqrt(M2_ST)*randn;
    LR_ST(i_reps,1)=0.5*(M1_ST - M2_ST + sqrt( (M1_ST- M2_ST)^2 + 4*M12_ST^2));
end

cv_LR = quantile(LR_ST, 1-alpha);