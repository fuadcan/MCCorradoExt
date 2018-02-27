output file=cza.out reset ;

   seed = 238475 ;
   dx =   rndns(100,2,seed) ;            /* Generate a random walk          */
   y  =   cumsumc(dx);                   /* Accumulate the innovations      */
                                         /* y is a non-cointegrated system  */
                                         /* of 2 I(1) variables             */

   _ker_fun = &parzen ;

   {alpha, t_1, t, c_a, c_t } = cza(y[.,1],y[.,2],-1,5) ;

   "True model is not cointegrated.";
   "PCB statistic for cointegration";
   "Autoregressive parameter = " alpha ;
   "PCB T*(alpha - 1); 5% critical value = " t_1~c_a[2,1] ;
   "PCB t-statistic; 5% critical value   = " t~c_t[2,1] ;
   "Accept null of no-cointegration."; " " ;

   /* Now develop a cointegrated system! */

   z = y[.,1] + dx[.,2] ;

   /* Z and Y are cointegrated.  Should easily reject null of no-coint  */

   "True model is cointegrated.";
   {alpha, t_1, t, c_a, c_t } = cza(z,y[.,1],-1,5) ;
   "PCB statistic for cointegration";
   "Autoregressive parameter = " alpha ;
   "PCB T*(alpha - 1); 5% critical value = " t_1~c_a[2,1] ;
   "PCB t-statistic; 5% critical value   = " t~c_t[2,1] ;
   "Reject null of no-cointegration." ;

