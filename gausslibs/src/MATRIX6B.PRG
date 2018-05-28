/*
** see arma1d.prg
*/

new;
library tsm,optmum;
TSMset;

/*
**  AR1_11 = AR1_21 (beta[1] = beta[2])
*/

CC = 1~-1~zeros(1,6);
c = 0;

{RR,r} = Explicit_to_Implicit(CC,c);
{CC2,c2} = Implicit_to_Explicit(RR,r);

format 6,3;
print "RR "; print RR;
print "r "; print r;
print "CC2 "; print CC2;
print "c2 "; print c2;

/*
** in arma1d.prg, we use the following instructions:
**
** w = 1|1|seqa(2,1,6);
** RR = design(w); r = zeros(8,1);
*/


