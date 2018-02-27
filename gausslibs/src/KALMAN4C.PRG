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
**                       _      _
**              _      _ | y(t) |
**     z(t)  =  | 1  0 | |      |  + e(t)
**              -      - | u(t) |
**                       -      -
**
**    _      _    _               _  _        _     _   _
**    | y(t) |    | phi1  -theta1 |  | y(t-1) |     | 1 |
**    |      |  = |               |  |        |  +  |   | u(t)
**    | u(t) |    |   0       0   |  | u(t-1) |     | 1 |
**    -      -    -               -  -        -     -   -
**
**  Maximum Likelihood in the frequency domain
*/

new;
library tsm,optmum,pgraph;

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


proc sgf_arma_plus_noise(beta,lambda);
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

sv = 0.95|0.5|sqrt(1|0.25);

_tsm_Mcov = 1;
_tsm_optmum = 1;
_print = 1;
_tsm_parnm = "phi1"|"theta1"|"sig_e"|"sig_u";
{theta1,stderr,Mcov,Logl} =  FD_ml(zt,&sgf_arma_plus_noise,sv);



theta2 = 0.9503|0.624|1.0828|0.000108;

{lambda,I} = PDGM(zt);
g1 = sgf_arma_plus_noise(theta1,lambda);
g2 = sgf_arma_plus_noise(theta2,lambda);

q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; _pnum = 2; fonts("simplex simgrma");
  title("ARMA(1,1) process with noise"\
        "\Lz]t[ = y]t[+e]t[   y]t[ = 0.95y]t-1[+u]t[-0.5u]t-1["\
        "\LVar@[e]t[@] = 0.25   Var@[u]t[@] = 2");
  _plwidth = 1|5|10;
  _plegstr = "Periodogram of z]t[\0s.g.f. with FDML parameters"\
             "\0s.g.f. with TDML parameters";
  _plegctl = {2 5 1 1};
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  xlabel("Frequency");
  asclabel(lab,0);
  logy(lambda[1:q],I[1:q]~g1[1:q]~g2[1:q]);







