output file=cadf.out reset ;

   seed = 238475 ;
   dx =   rndns(100,2,seed) ;            /* Generate a random walk          */
   y  =   cumsumc(dx) ;                  /* Accumulate the innovations      */
                                         /* y is a non-cointegrated system  */
                                         /* of 2 I(1) variables             */

   {alpha, t, c_t } = cadf(y[.,1],y[.,2],-1,5) ;

   "True model is not cointegrated.";
   "ADF statistic for cointegration";
   "Autoregressive parameter = " alpha ;
   "ADF t-statistic; 5% critical value   = " t~c_t[2,1] ;
   "Accept null of no-cointegration."; " " ;

   /* Now develop a cointegrated system! */

   z = y[.,1] + dx[.,2] ;

   /* Z and Y are cointegrated.  Should easily reject null of no-coint  */

   "True model is cointegrated.";
   {alpha,  t, c_t } = cadf(z,y[.,1],-1,5) ;
   "ADF statistic for cointegration";
   "Autoregressive parameter = " alpha ;
   "ADF t-statistic; 5% critical value   = " t~c_t[2,1] ;
   "Reject null of no-cointegration." ;

