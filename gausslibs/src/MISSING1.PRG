new;
library tsm,optmum;

rndseed 123;

x = rndu(20,1);
x[4] = miss(0,0);
x[12] = miss(0,0);

y = 2*x + rndn(20,1);
y[3] = miss(0,0);

{data,retcode} = MISSING(y~x,0);

y_ = data[.,1];
x_ = data[.,2];

beta = y_/x_;  /* OLS estimates */

yhat_ = x_*beta;
u_ = y_ - yhat_;

/* We use the retcode vector defined before */

{data,retcode} = MISSING(yhat_~u_,retcode);

yhat = data[.,1];
u = data[.,2];

print "              Obs            Y              Yhat          Residuals";
print;
data=seqa(1,1,20)~y~yhat~u;
call printfmt(data,1);







