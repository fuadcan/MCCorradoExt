output file=adf.out reset ;
library tsm, coint;
   seed = 238475 ;
   dx =   rndns(100,1,seed) ;              /* Generate a random walk          */
   y  =   cumsumc(dx) ;                    /* Accumulate the innovations      */
                                           /* Compute the ADF statistic       */

   {alpha,t, c_t } = adf(y,-1,5) ;
   "ADF statistic for random walk without drift";
   "Autoregressive parameter = " alpha ;
   "ADF t-statistic; 5% critical value   = " t~c_t[2,1] ;

   "Woops!; p < -1";
   {alpha, t, c_t } = adf(y,-2,5) ;

   "Woops!; p > 5; critical values not available";
   {alpha, t, c_t } = adf(y,6,5) ;
   "Autoregressive parameter = " alpha ;
   "ADF t-statistic; 5% critical value   = " t~c_t[2,1] ;



