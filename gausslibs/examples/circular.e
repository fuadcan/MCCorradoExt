
print "         circular.e  -  A circular regression model";
print;
print "   The log-likelihood for circular regression contains local";
print "   minima where the coefficients become very large.  When";
print "   that happens, restart from another point.";
print;


/* Simulate circular data */


  rndseed 34534567;
  nobs = 500;
  numInd = 2;

  mu = 0;
  k = 1;

  b = .5 * ones(numInd,1);
  x = rndu(nobs,numInd);
  y = rndvm(nobs,1,mu + gd(x*b),k);



  b0 = 5*ones(NumInd+2,1);

 output file = circular_out.txt reset;

 { b,f0,grd,ret } = QNewton(&lpr,b0);

  cov = invpd(hessp(&lpr,b));

  print "coefficients   standard errors";
  print;
  print b~sqrt(diag(cov));

 output off;


proc Gd(u);
   retp(2*atan2(2*u,2));  // atan2 is more accurate than atan across platforms
endp;

proc lpr(b);
    local dev;
    dev = y - b[2]- Gd(b[3] + x * b[3:rows(b)]);
    retp(sumc(ln(mbesselei0(exp(b[1]))) - exp(b[1]) * (cos(dev) - 1)));
endp;

print;
print;
print "Results saved in the file named circular_out.txt";
print "You will find this file in your working directory";
