output file=sw.out reset ;

   seed = 238475 ;
   dx =   rndns(100,2,seed) ;        /* Generate a random walk          */
   y  =   cumsumc(dx) ;              /* Accumulate the innovations      */
                                     /* y is a non-cointegrated system  */
                                     /* of 2 I(1) variables             */

   _ker_fun = &parzen ;

   "Non-cointegrated model.";
   {sw_stat, cvs} = sw(y,-1,5) ;

   "Stock-Watson statistics" ; sw_stat ;
   "Critical values = " cvs ;

   /* Now develop a cointegrated system! */

   z = y[.,1] + dx[.,2] ;

   /* Z and Y are cointegrated.  Should easily reject null of no-coint  */

   y = z~y[.,1] ;

   " " ; "Cointegrated model";

   {sw_stat, cvs} = sw(y,-1,5) ;

   "Stock-Watson statistics" ; sw_stat ;
   "Critical values = " cvs ;
