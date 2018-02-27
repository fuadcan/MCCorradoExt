new;
library tsm,pgraph;

cls;

output file = wt3.out reset;

N = 512;
x = rndn(N,1);

{H,G,Htilde,Gtilde} = Coiflet(6);

print "Coiflet Filters";
print "---------------";


i = 1;
do until i > 5;
  _wcenter = i;
  w = wt(x,H,G,0);
  print ftos(i,"delay = %lf",1,0)$+ftos(w[1]," ---  c0 value: %lf",10,5);
  i = i + 1;
endo;

print;
print "Pollen Filters";
print "---------------";

{H,G,Htilde,Gtilde} = Pollen(1.25,3.25);

i = 1;
do until i > 5;
  _wcenter = i;
  w = wt(x,H,G,0);
  print ftos(i,"delay = %lf",1,0)$+ftos(w[1]," ---  c0 value: %lf",10,5);
  i = i + 1;
endo;


print;
print "Daubechies Filters";
print "---------------";

{H,G,Htilde,Gtilde} = Daubechies(18);

i = 1;
do until i > 5;
  _wcenter = i;
  w = wt(x,H,G,0);
  print ftos(i,"delay = %lf",1,0)$+ftos(w[1]," ---  c0 value: %lf",10,5);
  i = i + 1;
endo;

M = ln(N)/ln(2);
c0 = sumc(x)/2^(M/2);
print; print ftos(c0,"sumc(x)/2^(M/2) = %lf",10,5);

output off;
