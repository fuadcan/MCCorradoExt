/*
** Monte Carlo investigations of the Hankel singular values
*/


new;
library tsm,optmum,pgraph;

svd_ = {};
Nobs = 100;
Ns = 1500;

i = 1;
do until i>Ns;

  y = rndn(Nobs,1);

  y = recserar(y,1|1,1.25|-0.5);
  H = Hankel(y,6);
  svd_ = svd_|svd(H)';

  i = i+1;
endo;

svd__ = svd_[.,2:cols(svd_)]./svd_[.,1];

freq = svd__ .<=  1.96/sqrt(Nobs);
freq = sumc(freq) / Ns;

print freq;


