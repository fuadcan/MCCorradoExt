cls;
print "Quantile Example";
print;

rndseed 345567;

x = rndn(1000,4);  /* data */

e = { .025, .5, .975 };  /* quantile levels */


y = quantile(x,e);


print "median " y[2,.];

print;

print "95 percentiles " y[1,.] "  " y[3,.];

