output file=jstat.out reset ;

   seed = 238475 ;
   dx =   rndns(100,1,seed) ;        /* Generate a random walk          */
   y  =   cumsumc(dx) ;
                                     /* Accumulate the innovations      */
                                     /* Compute the J(p,q) statistic    */

   {js, t_1} = jstat(y,-1,5) ;
   "J(p,q) statistic for random walk without drift";
   "J(-1,5); 5% critical value " js~t_1[2,1];
   "Accept null hypothesis of a unit root";

