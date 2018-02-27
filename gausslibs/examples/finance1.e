cls;


print "           GAUSS example finance1.e";
print;
print "This example computes the call premium and implied volatility";
print "for the OEX stock index option on January 30, 2001 with";
print "expiration date of February 16.  ";
print;
print "Press any key to continue"; call keyw;

     S0 = 718.46;
     K = { 720, 725, 730 };
     r = .0498;
     div = 0;
     sigma = .2493;
     t0 = dtday(2001,1,30);
     t1 = dtday(2001,2,16);
     tau = ElapsedTradingDays(t0,t1) / annualTradingDays(2001);

     c = EuropeanBSCall(S0,K,r,div,tau,sigma);
     print c; 
     c = { 13.70, 11.90, 9.10 };
     S0 = 718.46;
     K = { 720, 725, 730 };
     r = .0498;
     div = 0;
     t0 = dtday(2001,1,30);
     t1 = dtday(2001,2,16);
     tau = elapsedTradingDays(t0,t1) / annualTradingDays(2001);

     sigma = EuropeanBSCall_ImpVol(c,S0,K,r,div,tau);
     print sigma; 
     { d,g,t,v,rh } = EuropeanBSCall_Greeks(S0,K,r,div,tau,sigma);
     print "delta " d;
     print "gamma " g;
     print "theta " t;
     print "vega " v;
     print "rho " rh;

