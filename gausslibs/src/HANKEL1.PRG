/*
**  Hankel matrix of a univariate time series
*/


new;
library tsm,optmum;

rndseed 123456;

output file = hankel1.out reset;

Nobs = 1000;
y = recserar(rndn(Nobs,1),1,.95);
H = Hankel(y,4);

print "AR(1) process"; print chrs(45*ones(70,1));
print "Hankel Matrix: "; print H; print;
s = svd(H);
print "Singular values: "; print svd(H); print;
M = sumc(s./s[1].>1.96/sqrt(Nobs));
print ftos(M,"McMillan order: %lf",1,0);

print; print;

y = recserar(rndn(Nobs,1),0|0,.4|0.3);
H = Hankel(y,4);

print "AR(2) process"; print chrs(45*ones(70,1));
print "Hankel Matrix: "; print H; print;
s = svd(H);
print "Singular values: "; print svd(H); print;
M = sumc(s./s[1].>1.96/sqrt(Nobs));
print ftos(M,"McMillan order: %lf",1,0);

output off;


