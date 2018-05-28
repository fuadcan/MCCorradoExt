/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Autocovariances and autocorrelations of a VAR process
**
**  See LUTKEPOHL [1991], section 2.1.4, pages 21-27
*/

new;
library tsm,optmum;

let PHI1 = {.5 .1, .4 .5};
let PHI2 = {0 0, .25 0};
beta = vec(PHI1~PHI2);
let SIGMA = {.09 0, 0 .04};

call arma_autocov(beta,2,0,SIGMA,3);

/*
**  See LUTKEPOHL page 24 to compare the results
*/

i = 0;
do until i>3;
  nom_ = "acv"$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;

print; print; call pause(2);

/*
**  See LUTKEPOHL page 24 to compare the results
*/

i = 0;
do until i>3;
  nom_ = "acr"$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;





