output file=ps.out reset ;

   seed = 238475 ;
   d    = ones(100,1);
   dx   = rndns(100,2,seed) ;
   e    = rndns(100,1,seed) ;              /* Generate random walks           */
   x    = recserar(dx,0.0~0.0,1.0~1.0) ;   /* Accumulate the innovations      */
                                           /* x is a non-cointegrated system  */
                                           /* of 2 I(1) variables             */

   /* develop a cointegrated regression equation */
   /* cointegrating vector: (1,-1,-1)            */

   y  = 10.00 + x[.,1] + e ;

   _ker_fun = &parzen ;

   /* Estimate using ps: Pentti Saikkonen  */

   {beta, vc, stderr, sigma, tstats, rss, resid} = ps(y,x[.,1],d,4,4) ;

"True Regression Coefficients:";
 1.0|10.0 ;
" " ;
"P Saikkonen(1991) regression coefficients, tstats:";
 beta~tstats ;  " " ;


"Coefficient on constant (p = 0)          :" beta[1,1] ;
"Coefficient on x                         :" beta[2,1];
"Coefficient on lags(diff(x,1)) 1,2,3,4 :" beta[3:6,1];
"Coefficient on leads(diff(x,1)) 1,2,3,4:" beta[7:10,1];
"Coefficient on diff(x,1)               :" beta[11,1];


