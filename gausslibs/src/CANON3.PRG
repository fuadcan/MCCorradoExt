/*
**  LUTKEPOHL [1991], Introduction to Multiple Time Series Analysis,
**  Springer-Verlag, Berlin-Heidelberg
**
**  Moving Average Representation of a VAR process
**
**  See LUTKEPOHL [1991], pages 217-236
**
**  Compute the canonical form for the ARMA(1,1) process
**  in the example ARMA1.PRG
*/

new;
library tsm,optmum;

load beta = arma1;

call canonical_arma(beta,1,1,5);  /* This is a VARMA(1,1) process */

nom = "MA";
i = 0;
do until i>5;
  nom_ = nom$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;

print; print; call pause(2);

nom = "AR";
i = 0;
do until i>5;
  nom_ = nom$+ftos(i,"%lf",1,0);
  print $nom_;
  x=varget(nom_); print x;
  i = i+1;
endo;





