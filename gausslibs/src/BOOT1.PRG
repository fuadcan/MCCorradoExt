/*
** Bootstrapping a time-invariant state space model
** The SSM is a Vector ARMA model
** (see arma1a.prg)
*/

new;
library optmum,tsm,pgraph;
TSMset;

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

load theta = EXACT_ML;

beta = theta[1:8];
Pchol = (theta[9]~0)|(theta[10]~theta[11]);
SIGMA = Pchol*Pchol';

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);

{a0,P0} = SSM_ic;
call KFiltering(y,a0,P0);

{ys,as} = bootstrap_SSM(a0);

Nobs = rows(y);

graphset;
  title("Bootstrapping state space models (first series)");
  _pdate = ""; _pnum = 2;
  _plegstr = "Original data\0bootstrap data set";
  _plegctl = 1;
  s = seqa(1,1,Nobs);
  xy(s,y[.,1]~ys[.,1]);

  title("Bootstrapping state space models (second series)");
  xy(s,y[.,2]~ys[.,2]);



