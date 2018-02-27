new;
library tsm,optmum,pgraph;
TSMset;

phi1 = -0.9;
phi2 = -0.65;
theta1 = 0.80;
sigma = 2;

beta = phi1|phi2|theta1;

lambda = seqa(0,pi/50,51);
g1 = sgf_arfima(beta|0,2,1,sigma,lambda);

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,2,1,sigma^2);
call SSM_build(Z,d,H,T,c,R,Q,0);
g2 = sgf_SSM(lambda);


graphset;
  _pdate = ""; _pnum = 2;
  _ptitlht = 0.17; _pnumht = 0.15; _paxht = 0.15;
  fonts("simplex simgrma");
  title("spectral density function of an ARMA(2,1) process");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  _plegstr = "computed with sgf_arfima\0computed with sgf_arma";
  _plegctl = {2 5 1.5 4};
  xy(lambda,g1~g2);

