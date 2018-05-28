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


i = 0;
do until i>10;
  nom_ = "acv"$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;

print; print; call pause(2);

i = 0;
do until i>10;
  nom_ = "acr"$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;






