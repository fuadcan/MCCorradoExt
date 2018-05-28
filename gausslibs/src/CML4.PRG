/*
**  Exact Maximum Likelihood Estimation of an ARMA(1,1) model
**  with the Kalman Filter
*/

new;
library tsm,optmum,CML;

TSMset; CMLset;


/*  Generate an ARMA(1,1) process */

rndseed 123456;

et = rndn(101,1)*1.5;

yt = recserar(et[2:101]+0.6*et[1:100],0,0.8);


/*  Build a procedure for the likelihood function */

proc ml(coeff,data);
  local beta,sigma,Pchol,Z,d,H,T,c,R,Q;
  local a0,P0,Logl;

  beta = coeff[1:2];
  SIGMA = coeff[3]^2;

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  {a0,P0} = SSM_ic;
  call KFiltering(data,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

/*
**  coeff[1]>0.8
*/

proc ineq(coeff);
  retp(coeff[1]-0.8);
endp;

_CML_IneqProc = &ineq;

output file = CML4.out reset;

_CML_algorithm = 1;
{beta,f,g,cov,ret} = CMLPrt(CML(yt,0,&ml,0.9|-0.6|1.5));

output off;

