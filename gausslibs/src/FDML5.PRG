new;
library tsm,optmum,pgraph;
TSMset;

rndseed 123456;

_fourier = 0;

Z = eye(2); d = 0|0;
let H[2,2] = 0.2 0 0 0.1;
let T[2,2] = 0.5 0.3 0 -0.5;
c = 0|0; R = 1|1; Q = 0.25;

call SSM_build(Z,d,H,T,c,R,Q,0);
{y,a} = RND_SSM(0|0,128);

{lambda,Iy} = PDGM2(y);

sv = 0.5|0.3|-0.5|sqrt(0.25|0.2|0.1);

/*
** Procedure to compute the sgf of the SSM
*/

proc sgf(theta,lambda);
  local T,Q,H,Z,d,c,R;
  local g;

  T = (theta[1]~theta[2])|(0~theta[3]);
  Q = theta[4];
  Q = Q^2;
  H = (theta[5]~0)|(0~theta[6]);
  H = H^2;
  Z = eye(2); d = 0|0;
  c = {0,0}; R = {1,1};

  call SSM_build(Z,d,H,T,c,R,Q,0);
  g = sgf_SSM(lambda);
  retp(g);
endp;


proc ml(theta);
  local Gy,Nstar,logl,j,Ij,Gj;

  Gy = sgf(theta,lambda);
  Nstar = rows(lambda);

  logl = zeros(Nstar,1);

  j = 1;
  do until j > Nstar;
   Ij = xpnd2(Iy,j);
   Gj = xpnd2(Gy,j);
   logl[j] = -0.5*ln(det(Gj)) -0.5*sumc(diag((inv(Gj)*Ij)));
   j = j + 1;
  endo;

  logl = real(logl);

  retp(logl);
endp;

output file = fdml5.out reset;

_tsm_Mcov = 1;
{theta,stderr,Mcov,Logl} = TD_ml(&ml,sv);

output off;
