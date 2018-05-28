/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  This program reproduces the figures 3.8 and 3.9 of LUTKEPOHL [1991].
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
{theta,stderr,Mcov,LOGL} = varx_ML(data,1,2);  /* Constant and 2 lags */

/* Take the AR part of the VARX model */

beta = theta[1:18];
const = theta[19:21];  /* constant estimates */
SIGMA = _varx_SIGMA;

Nr = 8;        /* Number of responses */
Ns = 50;       /* Number of simulations */

/* Compute the orthogonal reponses */

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,0,SIGMA);
c[1:3] = const;   /* Add the constant */
call SSM_build(Z,d,H,T,c,R,Q,0);

{DELTA,KSI,ZETA} = SSM_orthogonal(0|1|0,Nr);
rep1 = DELTA[.,3];
rep2 = KSI[.,3];
rep3 = ZETA[3];

/* Compute the standard errors with the Bootstrap Method */

rep1s = zeros(1+Nr,Ns);
rep2s = zeros(1+Nr,Ns);

i = 1;
do until i > Ns;

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,0,SIGMA);
  c[1:3] = const;   /* Add the constant */
  call SSM_build(Z,d,H,T,c,R,Q,0);
  a0 = data[1,.]'|data[2,.]'|zeros(3,1);
  P0 = zeros(9,9);
  call KFiltering(data[3:Nobs,.],a0,P0);

  {ys,as} = bootstrap_SSM(a0);         /* SSM bootstrap */
  datas = data[1,.]|data[2,.]|ys;      /* Add the presample values */
  {thetas,stderr,Mcov,LOGL} = varx_LS(datas,1,2); /* VAR estimation */

  betas = thetas[1:18];
  CONSTs = thetas[19:21];
  SIGMAs = _varx_SIGMA;

  {Z,d,H,T,c,R,Q} = arma_to_SSM(betas,2,0,SIGMAs);
  c[1:3] = CONSTs;   /* Add the constant */
  call SSM_build(Z,d,H,T,c,R,Q,0);

  {DELTAs,KSIs,ZETAs} = SSM_orthogonal(0|1|0,Nr);

  j = 0;
  do until j > Nr;
    z =  DELTAs[1+j,.]';
    rep1s[1+j,i] = z[3];
    z =  KSIs[1+j,.]';
    rep2s[1+j,i] = z[3];
    j = j + 1;
  endo;

  i = i + 1;
endo;

stderr1 = stdc(rep1s');
stderr2 = stdc(rep2s');

period = seqa(0,1,Nr+1);

graphset;
  _pnum = 2; _pdate = ""; _pltype = 1|6|1; _pframe = {0,0}; _pcross = 1;
  _pcolor = 12|11|12;
  bound = rep1 + (-2~0~2).*stderr1;
  title("Fig. 3.8. Estimated responses of consumption to an orthogonalized"\
        "\Limpulse in income with two-standard error bounds");
  xy(period,bound);

graphset;
  _pnum = 2; _pdate = ""; _pltype = 1|6|1; _pframe = {0,0}; _pcross = 1;
  _pcolor = 12|11|12;
  bound = rep2 + (-2~0~2).*stderr2;
  _pline = 1~6~0~rep3~Nr~rep3~1~10~1;
  title("Fig. 3.9. Accumulated responses of consumption to an orthogonalized"\
        "\Limpulse in income with two-standard error bounds");
  xy(period,bound);


