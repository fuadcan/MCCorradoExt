/*
**  TONG [1990], Non-linear Time Series, Oxford University Press
**
**  PRIESTLEY [1994], Spectral Analysis and Times Series,
**  Academic Press, London
**
**  Lynx data (1821-1934)
**  Tong, page 470
*/

new;
library tsm,optmum,pgraph;

load lynx[] = lynx.asc;

y = log(lynx);
Nobs = rows(y);

t = seqa(1821,1,1934-1821+1);

graphset;                            /* Priestley, page 384 */
   title("Lynx data");
  _pdate = ""; _pnum = 2;
  xlabel("date");
  xy(t,y);


{lambda,I} = PDGM(y); I[1] = miss(0,0);

q = trunc(rows(lambda)/2);

graphset;                            /* Priestley, page 412 */
  _pdate = ""; _pnum = 2;
  fonts("simplex simgrma");
  title("Periodogram of the Lynx data");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xlabel("frequency");
  xy(lambda[1:q],I[1:q]);

/*
**  Fisher's g-statistic
**
**  Priestley, page 407, formula 6.1.66
*/

g = maxc(I[2:57])/sumc(I[2:57]);
gstar = 2*Nobs*g;                    /* Priestley, page 407, formula 6.1.64 */
pvalue = Nobs*exp(-gstar/2);         /* Priestley, page 407, formula 6.1.65 */

print "Test for Periodogram Ordinates";
print chrs(45*ones(60,1)); print;
s = maxindc(I[2:57]);
print ftos(lambda[s],"Test for periodic component at the frequency %lf",4,2);
print;
print ftos(g,     "Fisher's g-statistic : %lf",6,5);
print ftos(gstar, "g* :                   %lf",6,5);
print ftos(pvalue,"p-value :              %lf",6,5);

