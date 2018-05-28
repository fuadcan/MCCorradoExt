/*
**  Hankel matrix of multivariate time series
*/

new;
library tsm,optmum;
cls;
y = recserar(rndn(100,3),1~2~3,.95~.85~.62);
H = Hankel(y,4);
format 6,2;
print H;


