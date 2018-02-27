/*
**  HARVEY [1990], Forecasting, Structural Time Series and
**  the Kalman Filter, Cambridge University Press, pages 86-89
**
**  Cycle Model
*/

new;
library tsm,optmum,pgraph;

load rainfall[] = rainfall.asc;
Nobs = rows(rainfall);
t_ = seqa(1849,1,Nobs);

graphset;
  title("Rainfall");
  _pdate = "";
  xy(t_,rainfall);

y = rainfall[1:131];
y = y - 142.2;

_cycle_prmt = 0;
_tsm_parnm = "rho"|"lambda_c"|"sig_kappa"|"sig_epsilon";
_tsm_optmum = 1;
_print = 1;

/*
**  Using the starting values in Harvey, page 198
*/

rho = 0.5; lambda_c = 0.78; sig_epsilon = 1; sig_kappa = 1;
sv = rho|lambda_c|sig_epsilon|sig_kappa;

output file = cycle2.out reset;

{theta,stderr,Mcov,Logl} =  FD_ml(y,&sgf,sv);

print; print; call pause(2);

theta[3 4] = theta[3 4]^2;
beta = 0.84|0.41|200|1612;

print "  Parameter     Estimates       Harvey's Estimates";
print;
call printfmt(_tsm_parnm~theta~beta,0~1~111);
print; print "These are the same estimates because if we set:"; print;
print "       rho* = - rho     and     lambda* = -lambda_c-pi"; print;
print "Then we verify that the value taken by the log-likelihood "\
      "function is the same.";
_tsm_parnm = "rho*"|"lambda*"|"sig_kappa"|"sig_epsilon";
theta[1] = - theta[1];
theta[2] = -theta[2]-pi;
print;
print "  Parameter     Estimates       Harvey's Estimates";
print;
call printfmt(_tsm_parnm~theta~beta,0~1~111);

output off;


/*
**  Spectral generating function for a stochastic cycle plus noise model
*/

proc sgf(theta,lambda);
  local g_cycle,g_noise,g;
  g_cycle = _cycle_sgf(theta[1:3],lambda);
  g_noise = theta[4]^2;
  g = g_cycle + g_noise;
  retp(g);
endp;




