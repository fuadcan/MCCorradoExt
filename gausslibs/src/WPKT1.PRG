new;
library tsm;

cls;

rndseed 123;

output file = wpkt1.out reset;

format 15,5;

sigma = 5;  /* Try different values for sigma, for example 0.1 and 100 */
N = 16;
x = recserar(ones(N,1),10,0.5) + rndn(N,1)*sigma;

{H,G,Htilde,Gtilde} = Daubechies(14);

_wcenter = 7;
Cost = &Entropy;
{wpk1,B1} = BestBasis(x,H,G,Cost);
cost1 = _cost;
{wpk2,B2} = BestLevel(x,H,G,Cost);
cost2 = _cost;

print "       Best Basis   Best Level Basis";
m = maxc(rows(B1)|rows(B2));
B = miss(zeros(m,2),0);
B[1:rows(B1),1] = B1;
B[1:rows(B2),2] = B2;
print B;
print "Cost: ";
print cost1~cost2;
print;
print "Cost verification:";
print Entropy(wpk1)~Entropy(wpk2);

output off;

