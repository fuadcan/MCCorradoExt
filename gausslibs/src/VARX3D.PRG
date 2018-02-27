/*
** SIMULTANEOUS LINEAR MODELS
**
** JUDGE, G.G., R.C. HILL, W.E. GRIFFITHS, H. LUTKEPOHL and
** T.S. LEE [1988], Introduction to the Theory and Practice
** of Econometrics, John Wiley & Sons, New York
**
** See pages 636-663
*/

new;
library tsm,optmum;

cls;

let data[20,8] =

1        3.06    1.34    8.48    28      359.27  102.96  578.49
1        3.19    1.44    9.16    35      415.76  114.38  650.86
1        3.3     1.54    9.9     37      435.11  118.23  684.87
1        3.4     1.71    11.02   36      440.17  120.45  680.47
1        3.48    1.89    11.64   29      410.66  116.25  642.19
1        3.6     1.99    12.73   47      530.33  140.27  787.41
1        3.68    2.22    13.88   50      557.15  143.84  818.06
1        3.72    2.43    14.5    35      472.8   128.2   712.16
1        3.92    2.43    15.47   33      471.76  126.65  722.23
1        4.15    2.31    16.61   40      538.3   141.05  811.44
1        4.35    2.39    17.4    38      547.76  143.71  816.36
1        4.37    2.63    18.83   37      539     142.37  807.78
1        4.59    2.69    20.62   56      677.6   173.13  983.53
1        5.23    3.35    23.76   88      943.85  223.21  1292.99
1        6.04    5.81    26.52   62      893.42  198.64  1179.64
1        6.36    6.38    27.45   51      871     191.89  1134.78
1        7.04    6.14    30.28   29      793.93  181.27  1053.16
1        7.81    6.14    25.4    22      850.36  180.56  1085.91
1        8.09    6.19    28.84   38      967.42  208.24  1246.99
1        9.24    6.69    34.36   41      1102.61 235.43  1401.94;

X = data[.,1:5];
ENDOG = data[.,6:8];
EXOG = ENDOG~X;

w = 0|1|0|2|0|3|4|0|0|5|6|7|0|8|9|0|10|0|0|11|0|0|0|12;
RR = design(w);
r = zeros(24,1);


output file = varx3d.out reset;


/* OLS estimate, see JHGLL page 660 */

_print = 0;

SIGMA0 = zeros(3,3);
_varx_SIGMA = SIGMA0;

{theta,stderr,Mcov,Logl} = varx_CLS(ENDOG,EXOG,0,RR,r);

SIGMA1 = _varx_SIGMA;

G = theta[1:9];
G = reshape(theta,3,3) - eye(3);
B = theta[10:24]; B = reshape(B,5,3);
RF = -B*inv(G);  /* Reduced Form */
print; print; print chrs(43*ones(75,1));
print "OLS Estimates";
print chrs(43*ones(75,1));
print "GAMMA = "; print G;
print "B = "; print B;
print "PI = "; print RF;

print;

/* GLS-2S */

_print = 0;

_varx_SIGMA = SIGMA1;

{theta,stderr,Mcov,Logl} = varx_CLS(ENDOG,EXOG,0,RR,r);

SIGMA2 = _varx_SIGMA;

G = theta[1:9];
G = reshape(theta,3,3) - eye(3);
B = theta[10:24]; B = reshape(B,5,3);
RF = -B*inv(G);  /* Reduced Form */
print; print; print chrs(43*ones(75,1));
print "GLS-2S Estimates";
print chrs(43*ones(75,1));
print "GAMMA = "; print G;
print "B = "; print B;
print "PI = "; print RF;

print;

/* GLS-3S */

_print = 0;

_varx_SIGMA = SIGMA2;

{theta,stderr,Mcov,Logl} = varx_CLS(ENDOG,EXOG,0,RR,r);

SIGMA3 = _varx_SIGMA;

G = theta[1:9];
G = reshape(theta,3,3) - eye(3);
B = theta[10:24]; B = reshape(B,5,3);
RF = -B*inv(G);  /* Reduced Form */
print; print; print chrs(43*ones(75,1));
print "GLS-3S Estimates";
print chrs(43*ones(75,1));
print "GAMMA = "; print G;
print "B = "; print B;
print "PI = "; print RF;

print;

/* LIML */

_print = 0;

_varx_tol = 1e-20;

{theta,stderr,Mcov,Logl} = varx_CML(ENDOG,EXOG,0,RR,r);
G = theta[1:9];
G = reshape(theta,3,3) - eye(3);
B = theta[10:24]; B = reshape(B,5,3);
RF = -B*inv(G);  /* Reduced Form */
print; print; print chrs(43*ones(75,1));
print "LIML Estimates";
print chrs(43*ones(75,1));
print "GAMMA = "; print G;
print "B = "; print B;
print "PI = "; print RF;


output off;

