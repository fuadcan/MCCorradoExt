new;
library tsm,optmum,pgraph;
TSMset;

declare external r;

rndseed 123;
y = recserar(rndn(500,1),0|0,0.5|0.4);

_fourier = 1;
_tsm_optmum = 0;
r = 5;
sv = ones(r+1,1);

_tsm_parnm = "gamma1"|"gamma2"|"gamma3"|"gamma4"|"gamma5"|"sigma";

w = 1|2|3|4|0|5;
RR = design(w); r_ = zeros(6,1);
sv = invpd(RR'RR)*RR'*(sv-r_);

output file = fdml1b.out reset;

{theta,stderr,Mcov,Logl} = FD_cml(y,&sgf,sv,RR,r_);

G_matrix = vread(_ml_derivatives,"G_matrix");
I_matrix = vread(_ml_derivatives,"I_matrix");

LM = G_matrix'*inv(I_matrix)*G_matrix;
pvalue = cdfchic(LM,1);

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

{J_matrix,G_matrix,H_matrix,I_matrix} = FDml_derivatives(y,&sgf,theta);

LM = G_matrix'*inv(I_matrix)*G_matrix;
pvalue = cdfchic(LM,1);

print;
print ftos(LM,"Lagrange multiplier: %lf",10,5);
print ftos(pvalue,"p-value: %lf",10,5);

output off;


proc sgf(theta,lambda);       /* Dzhaparidze, page 125 */
  local N,gamma_,sigma,j,w,g;
  N = rows(lambda);
  gamma_ = theta[1:r];
  sigma = theta[r+1];
  j = seqa(1,1,r);
  w = lambda.*j';
  w = cos(w);
  g = (sigma^2)*exp(2*w*gamma_);
  retp(g);
endp;

