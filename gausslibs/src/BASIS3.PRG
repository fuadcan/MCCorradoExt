new;
library tsm;
TSMset;

x = rndn(8,1);

{H,G,Htilde,Gtilde} = Haar;

_wcenter = 0;
pkt = wpkt(x,H,G,0);

output file = basis3.out reset;

print "Packet table:"; print pkt; print;

{wpk,B} = BestLevel(x,H,G,&LogEnergy);
cost = _cost;

print "Basis: "; print B; print;
print "Wavelet Packet Coefficients:"; print wpk; print;
print ftos(cost,"Cost value: %lf",5,4);

B1 = 1|1;
B2 = 2|2|2|2;
B3 = 3|3|3|3|3|3|3|3;

wpk1 = Basis(pkt,B1);
cost1 = LogEnergy(wpk1);
print ftos(cost1,"Cost value associated with the B1 basis: %lf",5,4);

wpk2 = Basis(pkt,B2);
cost2 = LogEnergy(wpk2);
print ftos(cost2,"Cost value associated with the B2 basis: %lf",5,4);

wpk3 = Basis(pkt,B3);
cost3 = LogEnergy(wpk3);
print ftos(cost3,"Cost value associated with the B3 basis: %lf",5,4);



output off;
