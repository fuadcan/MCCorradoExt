output file=sj.out reset ;

   seed = 238475 ;
   dx =   rndns(100,2,seed) ;
   e  =   rndns(100,1,seed) ;              /* Generate random walks           */
   x  =   recserar(dx,0.0~0.0,1.0~1.0) ;   /* Accumulate the innovations      */
                                           /* x is a non-cointegrated system  */
                                           /* of 2 I(1) variables             */

   /* develop a cointegrated regression equation */
   /* cointegrating vector: (1,-1,-1)            */

   y  = 10.00 + x[.,1] + x[.,2] + e ;

   /* Estimate using SJ */
   {eval, evec, lr1, lr2} = sj(y~x,0,1) ;

"True Cointegrating vector:";
 1.00|-1.0|-1.0 ;
" " ;
" " ; "Johansen's MLE"; " " ;
"Eigen values:"; eval ;
" " ; "Cointegrating vector: "; evec[.,1] ;
" " ; "Johansen's trace statistics for 0 to" cols(y~x)-1 ":"  lr1 ;
" " ;
"Re-normalized cointegrating vector: " ; evec[.,1]/(evec[1,1]) ;
