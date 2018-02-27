output file=za.out reset ;

   _ker_fun = &parzen ;

   seed = 238475 ;
   dx   =   rndns(100,1,seed) ;        /* Generate a random walk       */
   y    =   cumsumc(dx) ;
                                     /* Accumulate the innovations      */
                                     /* Compute the PCB statistics      */

   {alpha, t_1, t, c_a, c_t } = za(y,-1,5) ;

   "PARZEN kernal ... " ;
   "PCB Z statistic for random walk without drift";
   "Autoregressive parameter = " alpha ;
   "PCB T*(alpha - 1); 5% critical value = " t_1~c_a[2,1] ;
   "PCB t-statistic; 5% critical value   = " t~c_t[2,1] ;

   "Woops!; p < -1";
   {alpha, t_1, t, c_a, c_t } = za(y,-2,5) ;

   "Woops!; p > 5; critical values not available";
   {alpha, t_1, t, c_a, c_t } = za(y,6,5) ;

   "Autoregressive parameter = " alpha ;
   "PCB T*(alpha - 1); 5% critical value = " t_1~c_a[2,1] ;
   "PCB t-statistic; 5% critical value   = " t~c_t[2,1] ;

    " " ; " " ; " " ; /* Now change to the FEJER kernal */
   _ker_fun = &fejer ;

   seed = 238475 ;
   dx =   rndns(101,1,seed) ;        /* Generate a random walk          */
   y  =   recserar(dx,0.0,1.0) ;     /* Accumulate the innovations      */
                                     /* Compute the ADF statistic       */

   {alpha, t_1, t, c_a, c_t } = za(y,-1,5) ;
   "FEJER kernal ... " ;
   "PCB Z statistic for random walk without drift";
   "Autoregressive parameter = " alpha ;
   "PCB T*(alpha - 1); 5% critical value = " t_1~c_a[2,1] ;
   "PCB t-statistic; 5% critical value   = " t~c_t[2,1] ;

   "Woops!; p < -1";
   {alpha, t_1, t, c_a, c_t } = za(y,-2,5) ;

   "Woops!; p > 5; critical values not available";
   {alpha, t_1, t, c_a, c_t } = za(y,6,5) ;
   "Autoregressive parameter = " alpha ;
   "PCB T*(alpha - 1); 5% critical value = " t_1~c_a[2,1] ;
   "PCB t-statistic; 5% critical value   = " t~c_t[2,1] ;

