/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Test for Instantaneous Causality
**  See LUTKEPOHL, section 3.6
**
*/

new;
library tsm,optmum;

cls;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUM = y[.,3];
dcons = CONSUM-lag1(CONSUM);
data = dinv~dinc~dcons;

output file = varx1c.out reset;

_print = 0;

/*
**  Use MLE to obtain the estimates of vech(chol(SIGMA))
*/

{theta,stderr,Mcov,LOGL} = varx_ML(data,1,2);  /* Constant and 2 lags */

/* theta is a 27*1 vector */

/*
**  Test for no Instantaneous-causality from income/consumption to investment:
**
**  beta[23] = beta[24] = 0
*/


RR = design(23|24)~zeros(2,3);

proc restrict(theta);
  retp(RR*theta);
endp;

/*
**  Note: LUTKEPOHL finds that the value of the statistic test is 5.70
*/

_print = 1;
{w,pvalue} = WALD(&restrict,theta,Mcov,rows(data));

output off;



/*
**  This procedure computes the Wald test
**
**  The hypothesis is: r(theta) = 0
**
**  r is the procedure r(theta)
**  T is the number of observations
*/

proc (2) = WALD(r,theta,Mcov,T);
  local r:proc;
  local r0,Dr0,g,w,pvalue;

  r0 = r(theta);
  g = rows(r0);       /* Number of restrictions */
  Dr0 = gradp(&r,theta);

  w = zeros(1,1); pvalue = w;

  /* Chi-2 test */

  w[1] = r0'*invpd(Dr0*Mcov*Dr0')*r0; pvalue[1] = cdfchic(w[1],g);


  if _print == 1;

    print ftos(w[1],"Wald test (chi-square): %lf",6,4);
    print ftos(pvalue[1],"p-value: %lf",6,4);

  endif;

  retp(w,pvalue);
endp;





