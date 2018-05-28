cls;
print;
print "This Example creates a random 600x600 positive definite matrix";
print "inverts it, and prints how long it took to complete these steps.";
print;

rndseed 34546;

output file = time_results.txt on;  

t0 = date;
x = rndu(700,600);
print "seconds to generate random matrix - " ethsec(t0,date)/100;

t0 = date;
y = moment(x,0);
print "seconds to compute cross-product - " ethsec(t0,date)/100;

t0 = date;
z = invpd(y);
print "seconds to compute inverse - " ethsec(t0,date)/100;

output off;

print;
print "Results appended to the time_results.txt file";
print "your working directory";
