new;
library tsm,optmum;

rndseed 123;

T = 1000;
_print = 1;
_fourier = 1; M = trunc(sqrt(T));
{y,retcode} = RND_arfima(0.5|0.5|0,1,1,2,1000,250,1);

/* Spectral estimation */

output file = gph5.out reset;

sv = 0.5|0.5|0|2;
cnt = 0|0|0|0;
call arfima(y,1,1,sv,cnt);

call pause(5);

/* Spectral estimation with d set equal to the GPH estimator */

{lambda,I} = PDGM(y);
d = GPH(I,lambda,M);
sv = 0.5|0.5|d|2;
cnt = 0|0|1|0;
call arfima(y,1,1,sv,cnt);

output off;

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



