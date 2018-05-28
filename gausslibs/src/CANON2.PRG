/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Moving Average Representation of a VAR process
**
**  See LUTKEPOHL [1991], pages 13-19
*/

new;
library tsm,optmum;

let PHI1 = {.5 .1, .4 .5};
let PHI2 = {0 0, .25 0};
beta = vec(PHI1~PHI2);

call canonical_arma(beta,2,0,3);  /* This is a VAR(2) process */

nom = "MA";
i = 0;
do until i>3;
  nom_ = nom$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;

print; print; call pause(2);

nom = "AR";
i = 0;
do until i>3;
  nom_ = nom$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;




