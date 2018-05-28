output file=fm-give.out reset ;

   seed = 238475 ;
   n    = 100 ;
   d    = ones(n,1);
   dx   = rndns(n,2,seed) ;
   dz   = rndns(n,4,seed) ;
   e2   = rndns(n,1,seed) ;
   e1   = rndns(n,1,seed) ;              /* Generate random walks           */
   x    = recserar(dx,0.0~0.0,1.0~1.0) ;   /* Accumulate the innovations      */
                                           /* x is a non-cointegrated system  */
                                           /* of 2 I(1) variables             */

   z    = recserar(dz,0.0~0.0~0.0~0.0,1.0~1.0~1.0~1.0) ;

   /* develop a cointegrated regression equation */
   /* cointegrating vector: (1,-1,-1)            */


   _times = 10 ;
   _filter = 0 ;
   _aband = 0;

   y1  = 0.2*x[.,1] + 0.5*x[.,2] + e1 ;
   y2  = 0.5*x[.,1] + 0.1*x[.,2] + e2 ;

   _ker_fun = &parzen;

   /* Estimate using FM-GMM */

   "system of two equations, two explanatory variables";
   {ahat, vc, s1,s2} = fm_give(y1~y2,x,z,3,3);
   "True coefficients = " 0.2~0.5 ;
   "                    " 0.5~0.1 ;
   "Ahat, FM-GIVE = " ahat ;
   "Long-Run Variance of the residuals " vc ;
   "s1 = " s1  ;
   "s2 = " s2  ;
   "st (= s1 + s2)" (s1+s2) ;
