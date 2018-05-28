/*
** see arma1e.prg
*/

new;
library tsm,optmum;
TSMset;

/*
**  AR1_12 = 0 (beta[3] = 0)
**  MA1_11 = 0 (beta[5] = 0)
**  MA1_21 = 0 (beta[6] = 0)
*/

CC = zeros(3,8);
CC[1,3] = 1;
CC[2,5] = 1;
CC[3,6] = 1;
c = 0|0|0;

{RR,r} = Explicit_to_Implicit(CC,c);

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

_print = 1;
_tsm_optmum = 0;
_tsm_gtol = 0.001;

output file = matrix6c.out reset;

{beta,stderr,Mcov,LogL} = arma_CML(y,1,1,0,RR,r);

output off;

/*
**  You can compare with the program arma1e.prg
**
**  Total observations:                                          100
**  Usable observations:                                          99
**  Number of parameters to be estimated:                          5
**  Degrees of freedom:                                           93
**  Value of the maximized log-likelihood function:       -503.30649
**
**
**  Parameters      estimates       std.err.      t-statistic         p-value
**  --------------------------------------------------------------------------
**  AR1_11          0.510686        0.087150         5.859863        0.000000
**  AR1_21          0.561029        0.125742         4.461759        0.000023
**  AR1_12          0.000000        0.000000                .               .
**  AR1_22          0.864802        0.040135        21.547077        0.000000
**  MA1_11          0.000000        0.000000                .               .
**  MA1_21          0.000000        0.000000                .               .
**  MA1_12         -0.158258        0.058780        -2.692366        0.008415
**  MA1_22          0.260918        0.106669         2.446062        0.016323
**
**  Covariance matrix: inverse of the negative Hessian.
**
*/


