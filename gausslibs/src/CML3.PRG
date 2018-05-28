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

output file = CML3.out reset;

{beta,f,g,cov,ret} = CMLPrt(CML(yt,0,&ml,0.8|-0.6|1.5));

output off;

/*
**
**  If you have the ARIMA library, you can compare:
**
**      {b,l,e,cv,aic,sbc} = arima(0,yt,1,0,1,0);
**
**  The results are:
**
**  Model:  ARIMA(1,0,1)
**
**
**  Final Results:
**
**  Iterations Until Convergence:   4
**
**  Log Likelihood:   -183.223266         Number of Residuals: 100
**  AIC           :    370.446532         Error Variance     : 2.294675936
**  SBC           :    375.656872         Standard Error     : 1.514818780
**
**  DF: 98      Adj. SSE: 228.550784245        SSE: 224.878241685
**
**              Coefficients     Std. Errors      T-Ratio    Approx. Prob.
**  AR1           0.71161142      0.07773737      9.15405          0.00000
**  MA1          -0.50727789      0.09713635     -5.22233          0.00000
**
**
**  Total Computation Time: 0.39 (seconds)
**
**  AR Root:   1.40526
**
**  MA Root:  -1.97131
*/

