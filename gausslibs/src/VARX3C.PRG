/*
** RESTRICTED SURE ESTIMATION
**
** JUDGE, G.G., R.C. HILL, W.E. GRIFFITHS, H. LUTKEPOHL and
** T.S. LEE [1988], Introduction to the Theory and Practice
** of Econometrics, John Wiley & Sons, New York
**
** See pages 460-462
*/

new;
library tsm,optmum;

cls;

let data[30,7] =

10.763   4.474   6.629   487.648   11.632  13.194    45.770
13.033  10.836  13.774   364.877   12.029   2.181    13.393
 9.244   5.856   4.063   514.037    8.196   5.586   104.819
 4.605  14.010   3.868   760.343   33.908   5.231   137.269
13.045  11.417  14.922   421.746    4.561  10.930    15.914
 7.706   8.755  14.138   578.214   17.594  11.854    23.667
 7.405   7.317   4.794   561.734   18.842  17.045    62.057
 7.519   6.360   3.768   301.470   11.637   2.682    52.262
 8.764   4.188   8.089   379.636    7.645  13.008    31.916
13.511   1.996   2.708   478.855    7.881  19.623   123.026
 4.943   7.268  12.901   433.741    9.614   6.534    26.255
 8.360   5.839  11.115   525.702    9.067   9.397    35.540
 5.721   5.160  11.220   513.067   14.070  13.188    32.487
 7.225   9.145   5.810   408.666   15.474   3.340    45.838
 6.617   5.034   5.516   192.061    3.041   4.716    26.867
14.219   5.926   3.707   462.621   14.096  17.141    43.325
 6.769   8.187  10.125   312.659    4.118   4.695    24.330
 7.769   7.193   2.471   400.848   10.489   7.639   107.017
 9.804  13.315   8.976   392.215    6.231   9.089    23.407
11.063   6.874  12.883   377.724    6.458  10.346    18.254
 6.535  15.533   4.115   343.552    8.736   3.901    54.895
11.063   4.477   4.962   301.599    5.158   4.350    45.360
 4.016   9.231   6.294   294.112   16.618   7.371    25.318
 4.759   5.907   8.298   365.032   11.342   6.507    32.852
 5.483   7.077   9.638   256.125    2.903   3.770    22.154
 7.890   9.942   7.122   184.798    3.138   1.360    20.575
 8.460   7.043   4.157   359.084   15.315   6.497    44.205
 6.195   4.142  10.040   629.378   22.240  10.963    44.443
 6.743   3.369  15.459   306.527   10.012  10.140    13.251
11.977   4.806   6.172   347.488    3.982   8.637    41.845;

data = ln(data);

p = data[.,1 2 3];   /* prices     */
y = data[.,4];       /* quantities */
q = data[.,5 6 7];   /* income     */

ENDOG = q;
EXOG = ones(30,1)~p~y;


output file = varx3c.out reset;

w = 1|2|3|4|0|0|0|5|0|0|0|6|7|8|9;
RR = design(w);
r = zeros(15,1);

/* First, perform a VARX estimation to obtain an estimate of SIGMA based
   on the unrestricted least squares residuals */

_print = 0;
_varx_sigma = 0;
call varx_cls(ENDOG,EXOG,0,RR,r);
SIGMA1 = _varx_sigma;

/* Or perform a VARX estimation to obtain an estimate of SIGMA based
   on the restricted least squares residuals */

w = 1|2|3|4|0|0|0|4|0|0|0|4|5|6|7;
RR = design(w);
r = zeros(15,1);

_varx_sigma = 0;
call varx_cls(ENDOG,EXOG,0,RR,r);
SIGMA2 = _varx_sigma;


/* Then, perform the Restricted SUR estimation */

_print = 1;

_varx_SIGMA = SIGMA1;
{theta1,stderr,Mcov,Logl} = varx_cls(ENDOG,EXOG,0,RR,r);

_varx_SIGMA = SIGMA2;
{theta2,stderr,Mcov,Logl} = varx_cls(ENDOG,EXOG,0,RR,r);


output off;






