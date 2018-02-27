output file=gstat.out reset ;

   seed = 238475 ;
   dx =   rndns(100,1,seed) ;        /* Generate a random walk          */
   y  =   cumsumc(dx) ;
                                     /* Accumulate the innovations      */
                                     /* Compute the G(p,q) statistic    */

   _ker_fun = &parzen ;

   {g, gp} = gstat(y,-1,5,10) ;
   "G(p,q) statistic for random walk without drift"; " ";
   "G(-1,5); p-value (PARZEN kernal)" g~gp;
   if ( gp < 0.05 ) ;
   "Null hypothesis of stationarity rejected using a 5% level of significance.";
   else ;
   "Null hypothesis accepted incorrectly";
   endif ;

" " ; "  " ;

   _ker_fun = &fejer ;

   {g, gp} = gstat(y,-1,5,10) ;
   "G(-1,5); p-value (FEJER kernal)" g~gp;
   if ( gp < 0.05 ) ;
   "Null hypothesis of stationarity rejected using a 5% level of significance.";
   else ;
     "Null hypothesis accepted --- incorrectly";
   endif ;

