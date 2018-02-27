output file=fm-ols.out reset ;

   seed = 238475 ;
   d    = ones(100,1);
   dx   = rndns(100,2,seed) ;
   e    = rndns(100,1,seed) ;              /* Generate random walks           */
   x    = recserar(dx,0.0~0.0,1.0~1.0) ;   /* Accumulate the innovations      */
                                           /* x is a non-cointegrated system  */
                                           /* of 2 I(1) variables             */

   /* develop a cointegrated regression equation */
   /* cointegrating vector: (1,-1,-1)            */



   _filter = 1 ;
   _aband = 0;
   _sbstart = 0.15 ;
   _sbend = 0.85 ;

   y  = 10.00 + x[.,1] + x[.,2] + e ;

   _ker_fun = &parzen ;

   /* Estimate using FM */

   _cregr   = &fm ;
   _hstat   = &hstatf ;

   {beta, vc} = fm_ols(y,x,d,5);

"True Regression Coefficients:";
 1.0|1.0|10.0 ;
" " ;
"Following estimates use PARZEN kernal"; " " ;
"FM-OLS regression coefficients, tstats:";
 beta ; vc;
beta./sqrt(diag(vc)) ;
output off ;
