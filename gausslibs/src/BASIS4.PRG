new;
library tsm;
declare external i;
TSMset;

x = recserar(rndn(4,1),0,1);

{H,G,Htilde,Gtilde} = Haar;

_wcenter = 0;
pkt = wpkt(x,H,G,0);

output file = basis4.out reset;

print "Packet table:"; print pkt; print;

_LpNorm = 0.5;

i = 1;
do until i > 3;

  print "---------------------------------------------";
  if i == 1;
    print "Best Basis with Entropy";
  elseif i == 2;
    print "Best Basis with LogEnergy";
  elseif i == 3;
    print "Best Basis with LpNorm";
  endif;
  print "---------------------------------------------";

  {wpk,B} = BestBasis(x,H,G,&cfunction);
  cost = _cost;

  print "Basis: "; print B; print;
  print "Wavelet Packet Coefficients:"; print wpk; print;
  print ftos(cost,"Cost value: %lf",5,4);

  B1 = 1|1;
  B2 = 1|2|2;
  B3 = 2|2|1;
  B4 = 2|2|2|2;

  wpk1 = Basis(pkt,B1);
  cost1 = cfunction(wpk1);
  print ftos(cost1,"Cost value associated with the B1 basis: %lf",5,4);

  wpk2 = Basis(pkt,B2);
  cost2 = cfunction(wpk2);
  print ftos(cost2,"Cost value associated with the B2 basis: %lf",5,4);

  wpk3 = Basis(pkt,B3);
  cost3 = cfunction(wpk3);
  print ftos(cost3,"Cost value associated with the B3 basis: %lf",5,4);

  wpk4 = Basis(pkt,B4);
  cost4 = cfunction(wpk4);
  print ftos(cost4,"Cost value associated with the B4 basis: %lf",5,4);


  i = i + 1;
endo;


output off;



proc cfunction(x);

  if i == 1;
    retp(Entropy(x));
  elseif i == 2;
    retp(LogEnergy(x));
  elseif i == 3;
    retp(LpNorm(x));
  endif;

endp;
