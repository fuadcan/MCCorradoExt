output file=ccr-fm.out reset ;

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

   {beta, vc, stderr, sigma, tstats, rss, resid, tests} = cregr(y,x,d,5);
   {h, hp} = hstatf(y,x,0,1,5) ;

"True Regression Coefficients:";
 1.0|1.0|10.0 ;
" " ;
"Following estimates use PARZEN kernal"; " " ;
"Phillips and Hansen(1990) FM regression coefficients, tstats:";
 beta~tstats ;  " " ;
"H(0,2) statistic, p-value : " h~hp ;

   _cregr   = &ccr ;
   _hstat   = &hstatc ;

   {beta, vc, stderr, sigma, tstats, rss, resid, tests} = cregr(y,x,d,5) ;
   {h, hp} = hstat(y,x,0,2,5) ;

" " ; " " ; "Park's (1992) CCR regression coefficients, tstats:";
 beta~tstats ;
" " ;
"H(0,2) statistic, p-value : " h~hp ;


" " ; " " ; "Now re-calculate using FEJER kernal"; " " ; " " ;

   /* Estimate using CCR */

   _ker_fun = &fejer ;
   _cregr   = &ccr ;
   _hstat   = &hstatc ;

   {beta, vc, stderr, sigma, tstats, rss, resid, tests} = cregr(y,x,d,5) ;
   {h, hp} = hstat(y,x,0,2,5) ;

"True Regression Coefficients:";
 1.0|1.0|10.0 ;  " " ;

"Park(1988) CCR regression coefficients, tstats:";
 beta~tstats ;  " " ;
"H(0,2) statistic, p-value : " h~hp ;

   _ker_fun = &fejer ;
   _cregr   = &fm ;
   _hstat   = &hstatf ;

   {beta, vc, stderr, sigma, tstats, rss, resid,tests} = cregr(y,x,d,5) ;
   {h, hp} = hstat(y,x,0,2,5) ;

" " ; " " ; "Phillips-Hansen FM regression coefficients, tstats:";
 beta~tstats ;
" " ;
"H(0,2) statistic, p-value : " h~hp ;
