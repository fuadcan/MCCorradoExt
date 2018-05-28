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

let PHI = {.5 0 0, .1 .1 .3, 0 .2 .3};
beta = vec(PHI);

call canonical_arma(beta,1,0,3);  /* This is a VAR(1) process */

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



