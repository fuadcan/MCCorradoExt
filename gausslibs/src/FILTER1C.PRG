/*
** Estimation of an univariate ARMA-GARCH process with
** the arma_Filter procedure
**
** The ARMA(2,1)-GARCH(1,1) model is
**
**     y(t) = phi1*y(t-1) + phi2*y(t_2) + epsilon(t) - theta1*epsilon(t-1)
**
**     epsilon(t) --> N(0,h(t)^2)
**
**     h(t)^2 = alpha0^2 + alpha1^2*epsilon(t-1)^2 + alpha2^2*h(t-1)^2
**
** In the example filter1b.prg, the hypothesis alpha2 = 0 can not be rejected.
** Suppose that alpha1 is 0 and that h(t) is constant. Then, we have
**
**         h(t)^2 = alpha0^2 + alpha2^2*h(t)^2
** or
**         h(t)^2 = sigma^2 = alpha0^2/(1-alpha2^2)
**
** In this example, we test the following hypothesis
**                alpha0 = sigma + alpha1 = 0 + alpha2 = 0
**
*/

new;
library tsm,optmum,pgraph;
TSMset;

Nobs = 200;
rndseed 123456;
y = RND_arma(0.6|0.3|0.25,2,1,1,miss(0,0),Nobs);

proc ml(beta);
  local e,h2,logL;

  e = arma_filter(y,beta[1:3],2,1);
  e = missrv(e,0);
  h2 = garch_filter(e,beta[4:6]^2,1,1);
  logL = -0.5*ln(2*pi)-0.5*ln(h2)-0.5*(e^2)./h2;

  retp(logL);
endp;

output file = filter1c.out reset;

_tsm_parnm  = "Phi1"|"Phi2"|"Theta1"|"alpha0"|"alpha1"|"alpha2";
_tsm_Mcov = 0;
sv = 0.6|0.3|0.25|1|0.5|0.5;
_print = 0;
{beta,stderr,Mcov,LogL} = TD_ml(&ml,sv);

sigma = sqrt(beta[4]^2/(1-beta[6]^2));
beta2 = beta[1:3]|sigma|0|0;

{J_matrix,G_matrix,H_matrix,I_matrix} = TDml_derivatives(&ml,beta2);

/*
**  Lagrange multiplier
**
**  DAVIDSON and MACKINNON [1990], Specification tests based on artificial
**  regressions, Journal of the American Statistical Association, 85, 220-227
**
**  Compute the LM test with an OPG artificial regression
*/

l = ones(Nobs,1);

J_matrix = J_matrix[.,1:4];
bols = l/J_matrix;                         /* OPG artificial regression */
resid = l-J_matrix*bols;                   /* Residuals                 */
R2 = 1 - (resid'resid)/(l'l);              /* R2 statistic              */
LM = Nobs*R2;                              /* T * R2                    */
pvalue = cdfchic(LM,3);

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;

/*
**  Note that you can obtain the same result with:
**
**  I = J_matrix'J_matrix;
**  LM = G_matrix'*inv(I)*G_matrix;
**
**  see for example:
**
**  DAVIDSON and MACKINNON [1993], Estimation and Inference in Econometrics
**  Oxford University Press, page 437, formula (13.04)
*/


