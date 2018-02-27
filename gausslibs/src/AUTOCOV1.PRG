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

let PHI = {0.5 0 0, .1 .1 .3, 0 .2 .3};
beta = vec(PHI);
let SIGMA = {2.25 0 0, 0 1.0 .5, 0 .5 .74};

call arma_autocov(beta,1,0,SIGMA,3);

/*
**  See LUTKEPOHL page 23 to compare the results
*/

i = 0;
do until i>3;
  nom_ = "acv"$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;

print; print; call pause(5);

/*
**  See LUTKEPOHL page 27 to compare the results
*/

i = 0;
do until i>3;
  nom_ = "acr"$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;




