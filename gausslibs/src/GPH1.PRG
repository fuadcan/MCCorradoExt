/*
**  GEWEKE and PORTER-HUDAK [1983], The estimation and application
**  of long memory time series model,
**  Journal of Time Series Analysis, 4, pages 221-238
*/

new;
library tsm,optmum;

cls;

T = 1000;               /* Number of observations              */

y = rndn(T,1);          /* Simulation of a white noise process */
_fourier = 0;           /* Using the fast Fourier transform    */
{lambda,I} = PDGM(y);

M = trunc(sqrt(T));     /* Using only the M lowest Fourier frequencies */

/*  See HURVICH and BELTRANO [1994], Journal of Time Series Analysis,
**  vol. 15 for a discussion of the choice of M
*/

Y = ln(I[2:M+1]);
w = 4*sin(lambda[2:M+1]/2)^2;
X = ones(M,1)~-ln(w);

beta = Y/X;
d = beta[2];

/* If you wish to compute the standard error of d */

XXinv = invpd(X'X);
std = pi/sqrt(6)*sqrt(XXinv[2,2]);

/* Testing d = 0 */

test = d/std; pvalue = 2*cdfnc(abs(test));

output file = gph1.out reset;

print ftos(d,"estimate of d     : %lf",5,4);
print ftos(std,"standard error    : %lf",5,4);
print ftos(test,"T-stat for d = 0  : %lf",5,4);
print ftos(pvalue,"p-value           : %lf",5,4);

output off;

