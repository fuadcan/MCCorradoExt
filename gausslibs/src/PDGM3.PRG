new;
library tsm,optmum,pgraph;

load sunspots[100,1] = sunspots.asc;

_fourier = 1;
{lambda,I} = PDGM(sunspots-meanc(sunspots));

format /rdn 16,4;
print I[1:50]';

print "Compare the results with these of Janacek and Swift, page 258";


q = trunc(rows(lambda)/2);

graphset;
  _pdate = ""; fonts("simplex simgrma"); _pnum = 2;
  title("Periodogram of the SUNSPOTS data");
  xlabel("frequency");
  xtics(0,pi,pi/4,0);
  lab = " 0 \202p\201/4 \202p\201/2 \2013\202p\201/4 \202p\201";
  asclabel(lab,0);
  xy(lambda[1:q],I[1:q]);

