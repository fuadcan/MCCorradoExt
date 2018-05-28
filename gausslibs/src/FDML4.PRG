/*
**  State space model of an ARMA(1,1) process with noise
**
**  Suppose that
**     z(t) = y(t) + e(t)
**  and
**     y(t) = phi1*y(t-1) + u(t) -theta1*u(t-1)
**
**  e(t) is the noise (the measure error for example)
**  y(t) is the ARMA(1,1) process
**  z(t) is the ARMA(1,1) process with noise
**
**  The state space form is:
**                       Ú      ¿
**              Ú      ¿ ³ y(t) ³
**     z(t)  =  ³ 1  0 ³ ³      ³  + e(t)
**              À      Ù ³ u(t) ³
**                       À      Ù
**
**    Ú      ¿    Ú               ¿  Ú        ¿    Ú   ¿
**    ³ y(t) ³    ³ phi1  -theta1 ³  ³ y(t-1) ³    ³ 1 ³
**    ³      ³  = ³               ³  ³        ³ +  ³   ³ u(t)
**    ³ u(t) ³    ³   0       0   ³  ³ u(t-1) ³    ³ 1 ³
**    À      Ù    À               Ù  À        Ù    À   Ù
**
**  Maximum Likelihood in the frequency domain
*/

new;
library tsm,optmum,pgraph;
TSMset;

/*
**  Simulation of the process with
**
**  var[e(t)] = 0.25, phi1 = 0.95, theta1 = 0.5 and var[u(t)] = 1
**  y(0) = 40;
*/

rndseed 123456;

t_ = seqa(1,1,500);

u = rndn(500,1)*sqrt(1);     /* Simulate the u(t) process */
u_ = u~(0|trimr(u,0,1));
u_ = u_*(1|-0.5);
y = recserar(u_,40,0.95);   /* Simulate the y(t) process */
e = rndn(500,1)*sqrt(0.25);
zt = y + e;                 /* Simulate the z(t) process */


proc sgf1(beta,lambda);
  local phi1,theta1,sig_u,sig_e;
  local w,w1,w2,g;
  phi1 = beta[1];
  theta1 = beta[2];
  sig_u = beta[3];
  sig_e = beta[4];
  w = cos(lambda);
  w1 = 1 - 2*theta1*w + theta1^2;
  w2 = 1 - 2*phi1*w + phi1^2;
  g = (w1./w2)*sig_u^2 + sig_e^2;
  retp(g);
endp;


/*
** Procedure to compute the sgf of the SSM
*/

proc sgf2(beta,lambda);
  local phi1,theta1,sig_u,sig_e;
  local T,Q,H,Z,d,c,R;
  local G;

  phi1 = beta[1];
  theta1 = beta[2];
  sig_u = beta[3];
  sig_e = beta[4];
  T = (phi1~-theta1)|(0~0);
  Q = sig_u^2; H = sig_e^2;
  Z = {1 0}; d = 0;
  c = {0,0}; R = {1,1};

  call SSM_build(Z,d,H,T,c,R,Q,0);
  G = sgf_SSM(lambda);
  G = real(G);
  retp(G);
endp;

sv = 0.95|0.5|sqrt(1|0.25);

_tsm_Mcov = 0;
_tsm_optmum = 1;
_tsm_parnm = "phi1"|"theta1"|"sig_e"|"sig_u";

{theta1,stderr,Mcov,Logl} =  FD_ml(zt,&sgf1,sv);

{theta2,stderr,Mcov,Logl} =  FD_ml(zt,&sgf2,sv);

print "             FDML with sgf1   FDML with sgf2";
call printfmt(_tsm_parnm~theta1~theta2,0~1~1);
print;
print "The first estimation is much faster, because the second sgf is
computed within a loop.";
