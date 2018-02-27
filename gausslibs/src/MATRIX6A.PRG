/*
** see arma1c.prg
*/

new;
library tsm,optmum;
TSMset;

/*
**  AR1_11 = 1 (beta[1] = 1)
*/

CC = 1~zeros(1,7);
c = 1;

{RR,r} = Explicit_to_Implicit(CC,c);
{CC2,c2} = Implicit_to_Explicit(RR,r);

format 6,3;
print "RR "; print RR;
print "r "; print r;
print "CC2 "; print CC2;
print "c2 "; print c2;

/*
** in arma1c.prg, we use the following instructions:
**
** RR = design(0|seqa(1,1,7)); r = 1|zeros(7,1);
*/

