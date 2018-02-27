/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  This program reproduces the table 3.3 of LUTKEPOHL [1991].
*/

new;
library tsm,optmum,pgraph;
TSMset;

load y[92,3] = lutkepoh.asc;

y = ln(y[1:76,.]);

INVEST = y[.,1];
dinv = INVEST-lag1(INVEST);
INCOME = y[.,2];
dinc = INCOME-lag1(INCOME);
CONSUM = y[.,3];
dcons = CONSUM-lag1(CONSUM);
{data,retcode} = Missing(dinv~dinc~dcons,0);
Nobs = rows(data);

/* First, estimate the coefficients */

_print = 0;
{theta,stderr,Mcov,LOGL} = varx_LS(data,1,2);  /* Constant and 2 lags */

/* Take the AR part of the VARX model */

beta = theta[1:18];
const = theta[19:21];  /* constant estimates */
SIGMA = _varx_SIGMA;

Nr = 8;        /* Number of responses */
Ns = 50;       /* Number of simulations */

rep1  = zeros(1+Nr,9);  /* Estimated responses */
rep2  = zeros(1+Nr,9);  /* Accumulated responses */

/* Compute the responses */

call arma_impulse(beta,2,0,Nr);

j = 0;
do until j > Nr;
  z = varget("IMPULSE"$+ftos(j,"%lf",1,0));
  rep1[1+j,.] = vec(z)';
  z = varget("_IMPULSE"$+ftos(j,"%lf",1,0));
  rep2[1+j,.] = vec(z)';
  j = j + 1;
endo;

/* Build the corresponding state space model */

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,0,SIGMA);
c[1:3] = const;   /* Add the constant */
call SSM_build(Z,d,H,T,c,R,Q,0);
a0 = data[1,.]'|data[2,.]'|zeros(3,1);
P0 = zeros(9,9);
call KFiltering(data[3:Nobs,.],a0,P0);

/* Compute the standard errors with the Bootstrap Method */

stderr1 = zeros(1+Nr,9);
stderr2 = zeros(1+Nr,9);

i = 1;
do until i > Ns;

  {ys,as} = bootstrap_SSM(a0);         /* SSM bootstrap */
  datas = data[1,.]|data[2,.]|ys;      /* Add the presample values */
  {thetas,stderr,Mcov,LOGL} = varx_LS(datas,1,2); /* VAR estimation */
  betas = thetas[1:18];
  call arma_impulse(betas,2,0,Nr);

  j = 0;
  do until j > Nr;
    z =  varget("IMPULSE"$+ftos(j,"%lf",1,0));
    stderr1[1+j,.] = stderr1[1+j,.] + (vec(z)'-rep1[1+j,.])^2;
    z =  varget("_IMPULSE"$+ftos(j,"%lf",1,0));
    stderr2[1+j,.] = stderr2[1+j,.] + (vec(z)'-rep2[1+j,.])^2;
    j = j + 1;
  endo;

  i = i + 1;
endo;

stderr1 = sqrt(stderr1/Ns);
stderr2 = sqrt(stderr2/Ns);

period = seqa(0,1,Nr+1);

let fmt[3,3] =
               "*.*lf" 10 3
               "*.*lf" 10 3
               "*.*lf" 10 3;

output file = impuls1d.out reset;

i = 0;
do until i > Nr;

  x = xpnd2(rep1,1+i);
  y = xpnd2(stderr1,1+i);

  print chrs(45*ones(50,1));
  print ftos(i,"Periods: %lf",1,0);
  print chrs(45*ones(50,1));
  print "        Estimated responses    ";
  print "           THETA matrix        ";
  print chrs(45*ones(50,1));
  call printfm(x,1,fmt);
  print;
  print "     Estimated standard errors  ";
  print chrs(45*ones(50,1));
  call printfm(y,1,fmt);
  print;

  i = i+1;
endo;

i = 0;
do until i > Nr;

  x = xpnd2(rep2,1+i);
  y = xpnd2(stderr2,1+i);

  print chrs(45*ones(50,1));
  print ftos(i,"Periods: %lf",1,0);
  print chrs(45*ones(50,1));
  print "        Accumulated responses    ";
  print "             PSI matrix        ";
  print chrs(45*ones(50,1));
  call printfm(x,1,fmt);
  print;
  print "     Estimated standard errors  ";
  print chrs(45*ones(50,1));
  call printfm(y,1,fmt);
  print;

  i = i+1;
endo;

output off;

