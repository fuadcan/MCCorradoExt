output file=pcb.out reset ;

   seed = 238475 ;
   d  =   ones(100,1) ;
   dx =   rndns(100,2,seed) ;
   e  =   rndns(100,1,seed) ;              /* Generate random walks           */
   x  =   recserar(dx,0.0~0.0,1.0~1.0) ;   /* Accumulate the innovations      */
                                           /* x is a non-cointegrated system  */
                                           /* of 2 I(1) variables             */

   /* develop a cointegrated regression equation */
   /* cointegrating vector: (1,-1,-1)            */

   y  = 10.00 + x[.,1] + x[.,2] + e ;

   /* Estimate using PCB */

{beta, vc, stderr, tstats} = pcb(y,x,d,2) ;
"True Regression Coefficients:";
 1.0|1.0|10.0;
" " ;
"Phillips spectral regression coefficients using PCB:";
 beta~tstats ;  " " ;

   /* Now estimate using B(0) */

{beta, vc, stderr, tstats} = pcbz(y,x,d,2) ;

" " ;
"Phillips spectral regression coefficients using PCBZ:";
 beta~tstats ;  " " ;
