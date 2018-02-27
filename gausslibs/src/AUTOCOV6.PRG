new;
library tsm,optmum,pgraph;
TSMset;

/* Time-invariant SSM with y(t) = [ y1(t) y2(t) ]' */

Z = eye(2); d = 10|0;
let H[2,2] = 2 0 0 1;
let T[2,2] = 0.5 0.3 0 0.2;
c = 1|4; R = 1|1; Q = 1;

call SSM_build(Z,d,H,T,c,R,Q,0);   /* cn = 0 */

{cv,cr} = SSM_autocov(10);
lag_ = seqa(0,1,11);

/*
** cv[.,1] ---> E[y1(t)y1(t-lag_)]
** cv[.,2] ---> E[y2(t)y1(t-lag_)]
** cv[.,3] ---> E[y1(t)y2(t-lag_)]
** cv[.,4] ---> E[y2(t)y2(t-lag_)]
*/

graphset;
  begwind;
  makewind(9,1,0,5.855,0);
  makewind(9/2,5.855/2,0,5.855/2,0);
  makewind(9/2,5.855/2,4.5,5.855/2,0);
  makewind(9/2,5.855/2,0,0,0);
  makewind(9/2,5.855/2,4.5,0,0);

  setwind(1);
    title("Autocorrelations of the State Space Model");
    _paxes = 0;
    _pnum = 0;
    _ptitlht = 1;
   draw;
  setwind(2);
    graphset;
    _ptitlht = 0.35; _paxht = 0.25; _pnumht = 0.25;
    _pnum = 2;
    title("E@[y]1[(t)y]1[(t-lag_)@]");
    xtics(-1,10,1,0);
    _pbarwid = 0.25; _pbartyp = {2 9};
    bar(lag_,cr[.,1]);
  setwind(3);
    title("E@[y]2[(t)y]1[(t-lag_)@]");
    bar(lag_,cr[.,2]);
  setwind(4);
    title("E@[y]1[(t)y]2[(t-lag_)@]");
    bar(lag_,cr[.,3]);
  setwind(5);
    title("E@[y]2[(t)y]2[(t-lag_)@]");
    bar(lag_,cr[.,4]);
  endwind;




