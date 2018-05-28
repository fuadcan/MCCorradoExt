/*
** see arma1f.prg
*/

new;
library tsm,optmum;
TSMset;

/*
**  AR1_21 = 0 (beta[2] = 0)
**  AR1_12 = 0 (beta[3] = 0)
**  MA1_21 = 0 (beta[6] = 0)
**  MA1_12 = 0 (beta[7] = 0)
*/

CC = zeros(4,8);
CC[1,2] = 1;
CC[2,3] = 1;
CC[3,6] = 1;
CC[4,7] = 1;

c = 0|0|0|0;

{RR,r} = Explicit_to_Implicit(CC,c);

load reinsel[100,2] = reinsel.asc;

invest = reinsel[.,1];
invent = reinsel[.,2];
di = invest - lag1(invest);
y = di~invent;

_print = 1;
_tsm_optmum = 0;
_tsm_gtol = 0.001;

output file = matrix6d.out reset;

{beta,stderr,Mcov,LogL} = arma_CML(y,1,1,0,RR,r);

output off;

/*
**  You can compare with the program arma1f.prg
**
**  Total observations:                                          100
**  Usable observations:                                          99
**  Number of parameters to be estimated:                          4
**  Degrees of freedom:                                           94
**  Value of the maximized log-likelihood function:       -513.38787
**
**
**  Parameters      estimates       std.err.      t-statistic         p-value
**  --------------------------------------------------------------------------
**  AR1_11          0.284297        0.182115         1.561083        0.121865
**  AR1_21          0.000000        0.000000                .               .
**  AR1_12          0.000000        0.000000                .               .
**  AR1_22          0.895241        0.051069        17.529942        0.000000
**  MA1_11         -0.266600        0.187651        -1.420724        0.158706
**  MA1_21          0.000000        0.000000                .               .
**  MA1_12          0.000000        0.000000                .               .
**  MA1_22          0.198175        0.114039         1.737787        0.085523
**
**  Covariance matrix: inverse of the negative Hessian.
**
*/



