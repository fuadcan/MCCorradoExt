/*
**  For a univariate ARMA process, we can use canonical_arma or
**  canonical_arfima
**
**  Suppose that we have:
**
**     y(t)-0.5*y(t-1)-0.25*y(t-2) = u(t)-0.4*u(t-1)+0.3*u(t-2)
*/

new;
library tsm,optmum;

beta = 0.5|0.25|0.4|-0.3|0;
{phi1,theta1,retcode} = canonical_arfima(beta,2,2,0,6);

call canonical_arma(beta[1:4],2,2,6);

phi2 = {}; theta2 = {};

nom = "MA";
i = 0;
do until i>6;
  nom_ = nom$+ftos(i,"%lf",1,0);
  x=varget(nom_); theta2 = theta2|x;
  i = i+1;
endo;

nom = "AR";
i = 0;
do until i>6;
  nom_ = nom$+ftos(i,"%lf",1,0);
  x=varget(nom_); phi2 = phi2|x;
  i = i+1;
endo;

print "AR representation"; print chrs(45*ones(40,1));
print  "   canonical_arfima    canonical_arma";
call printfmt(phi1~phi2,1~1);
print;
print "MA representation"; print chrs(45*ones(40,1));
print  "   canonical_arfima    canonical_arma";
call printfmt(theta1~theta2,1~1);




