new;
library tsm,optmum;

cls;

format /m1 /rd 2,0;
L = elimination_(4);
print "Elimination matrix L(4) : ";
print l;

D = duplication_(4);
print;
print "Duplication matrix D(4) : ";
print D;

K = commutation_(4,3);
print;
print "Commutation matrix K(4,3) : ";
print K;

