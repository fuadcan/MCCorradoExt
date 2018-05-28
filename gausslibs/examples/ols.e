
rndseed 345346;
x = rndn(100,5);
y = rndu(100,1);

output file = ols_output.txt reset;

call ols("",y,x);

output off;

print;
print "Results of this run stored in ols_output.txt file";
print "The file will be found in your working directory";

