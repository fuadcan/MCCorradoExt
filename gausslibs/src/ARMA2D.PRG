/*
**  Exact Maximum Likelihood Estimation of an ARMA(1,1) model
**  with the Kalman Filter
*/

new;
library tsm,optmum;
TSMset;

/*  Generate an ARMA(1,1) process */

rndseed 123456;

Nobs = 100;
et = rndn(Nobs+1,1)*1.5;
yt = recserar(et[2:Nobs+1]+0.6*et[1:Nobs],0,0.8);


/*  Build a procedure for the likelihood function */

proc ml(coeff);
  local beta,sigma,Pchol,Z,d,H,T,c,R,Q;
  local a0,P0,Logl;

  beta = coeff[1:2];
  SIGMA = coeff[3]^2;

  {Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
  call SSM_build(Z,d,H,T,c,R,Q,0);
  {a0,P0} = SSM_ic;
  call KFiltering(yt,a0,P0);
  Logl = KF_ml;

  retp(Logl);
endp;

_print = 1;
_tsm_optmum = 0;
_tsm_gtol = 0.001;
_tsm_Mcov = 1;

output file = arma2d.out reset;

{coeff,stderr,Mcov,Logl} = TD_ml(&ml,0.8|-0.6|1.5);

beta = coeff[1:2];
SIGMA = coeff[3]^2;

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);
{a0,P0} = SSM_ic;
call KFiltering(yt,a0,P0);
np = 10;
{Forecasts,mse,aF,PF} = KForecasting(np);
s = seqa(Nobs+1,1,np);
Fstderr = sqrt(mse);
LCL = Forecasts - 1.96*Fstderr;
UCL = Forecasts + 1.96*Fstderr;

omat = LCL~Forecasts~UCL~Fstderr;

print;
print "   Period        LCL         Forecasts        UCL"\
      "        Forecast Std. Err.";
print chrs(45*ones(75,1));
let fmt[5,3] = "*.*lf" 8 0
               "*.*lf" 15 6
               "*.*lf" 15 6
               "*.*lf" 15 6
               "*.*lf" 15 6;
call printfm(s~omat,1~1~1~1~1,fmt);

output off;

/*
**
**  If you have the ARIMA library, you can compare:
**
**      {b,l,e,cv,aic,sbc} = arima(0,yt,1,0,1,0);
**      f = forecast(b,yt,1,0,1,0,e,10);
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
**  Total Computation Time: 0.11 (seconds)
**
**  AR Root:   1.40526
**
**  MA Root:  -1.97131
**
**
**  Forecasts for ARIMA(1,0,1) Model.   95% Confidence Interval Computed.
**
**   Period        LCL        Forecasts        UCL     Forecast Std. Err.
**     101     -2.070068      0.898923      3.867913      1.514819
**     102     -4.041248      0.639684      5.320616      2.388275
**     103     -4.887352      0.455206      5.797765      2.725845
**     104     -5.324186      0.323930      5.972046      2.881745
**     105     -5.566196      0.230512      6.027220      2.957559
**     106     -5.706484      0.164035      6.034554      2.995218
**     107     -5.790815      0.116729      6.024274      3.014109
**     108     -5.843140      0.083066      6.009272      3.023630
**     109     -5.876523      0.059111      5.994744      3.028440
**     110     -5.898338      0.042064      5.982465      3.030873
*/
