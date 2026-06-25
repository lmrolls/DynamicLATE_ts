function b=ols_simple(y,x)

b=inv(x'*x)*x'*y;
