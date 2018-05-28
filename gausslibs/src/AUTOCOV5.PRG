/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Autocovariances and autocorrelations of a VARMA(1,1) process
**
**  See LUTKEPOHL [1991], section 6.4, pages 226-228
*/

new;
library tsm,optmum;

load beta = arma1;
load SIGMA = sigma1;

call arma_autocov(beta,1,1,SIGMA,10);  /* Vector ARMA(1,1) */

{Z,d,H,T,c,R,Q} = arma_to_SSM(beta,1,1,SIGMA);
call SSM_build(Z,d,H,T,c,R,Q,0);
{cv,cr} = SSM_autocov(10);

output file = autocov5.out reset;

i = 0;
do until i > 10;
  nom_ = "acv"$+ftos(i,"%lf",1,0);
  print; print $nom_;
  print "      computed with arma_autocov         computed with sgf_SSM    "  ;
  x = varget(nom_);
  Gi = xpnd2(cv,i+1);

  print x~Gi;

  i = i + 1;
endo;

i = 0;
do until i > 10;
  nom_ = "acr"$+ftos(i,"%lf",1,0);
  print; print $nom_;
  print "      computed with arma_autocov         computed with sgf_SSM    "  ;
  x = varget(nom_);
  Ri = xpnd2(cr,i+1);

  print x~Ri;

  i = i + 1;
endo;

output off;

