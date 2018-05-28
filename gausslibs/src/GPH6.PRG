/*
**  Monte Carlo investigations of the the GPH estimator's standard error
*/

new;
library tsm,optmum,pgraph;

T = 1000;
Ns = 250; /* Number of simulations */

_fourier = 0; M = trunc(sqrt(T)); d = zeros(Ns,3);
_smoothing = 2|0|1|0.20; /* a = .20 */

j = 1;
do until j>Ns;

  y = rndn(T,1);
  {lambda,I} = PDGM(y);

  d1 = GPH(I,lambda,M);      /* Geweke-Porter-Hudak     */

  _smoothing[2] = 2;         /* Bartlett                */
  Is1 = smoothing(I);
  d2 = GPH(Is1,lambda,M);

  _smoothing[2] = 4;         /* Tukey with a = .20      */
  Is2 = smoothing(I);
  d3 = GPH(Is2,lambda,M);

  d[j,.] = d1~d2~d3;

  j = j+1;
endo;



/*
** Procedure for the GPH estimator
*/

proc (1) = gph(I,lambda,M);
  local Y,w,X,beta,d;

  Y = ln(I[2:M+1]);
  w = 4*sin(lambda[2:M+1]/2)^2;
  X = ones(M,1)~-ln(w);

  beta = Y/X;
  d = beta[2];

  retp(d);
endp;


/* Variances of the smoothed periodogram */

var = zeros(3,1);
var[1] = pi^2/6;                /* periodogram              */
var[2] = 2/3;                   /* Bartlett                 */
var[3] = 2*(1-4*0.20+6*0.20^2); /* Tukey with a = 0.20;     */

w = 4*sin(lambda[2:M+1]/2)^2;
X = ones(M,1)~-ln(w);
XXinv = invpd(X'X);
std_ = sqrt(XXinv[2,2]*var);

output file = gph6.out reset;

let names = GPH Bartlett Tukey;
print "Standard error estimates"; print chrs(45*ones(40,1));
print "              Monte Carlo ";
call printfmt(names~stdc(d)~std_,0~1~1);

output off;
