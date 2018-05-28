output file=pupz.out reset ;

   seed = 238475 ;
   dx =   rndns(100,2,seed) ;            /* Generate a random walk          */
   y  =   cumsumc(dx) ;                  /* Accumulate the innovations      */
                                         /* y is a non-cointegrated system  */
                                         /* of 2 I(1) variables             */

   _ker_fun = &parzen ;

   {pu_stat, pu_crit} = pu(y[.,1],y[.,2],-1,5) ;

   "True model is not cointegrated.";
   "PU statistic for no-cointegration" pu_stat ;
   "5% critical value   = " pu_crit[5,1] ;
   "Accept null of no-cointegration."; " ";

   /* Now develop a cointegrated system! */

   z = y[.,1] + dx[.,2] ;

   /* Z and Y are cointegrated.  Should easily reject null of no-coint  */

   {pz_stat, pz_crit} = pz(z,y[.,1],-1,5) ;

   "True model is co-integrated.";
   "PZ statistic for no-cointegration" pz_stat ;
   "5% critical value   = " pz_crit[5,1] ;
   "Reject null of no-cointegration." ;
