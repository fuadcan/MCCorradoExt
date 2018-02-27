/*
**  The Kolmogorov-Smirnov test for white noise processes
**
**  BROCKWELL and DAVIS [1991], Times Series: Theory and Methods,
**  Springer-Verlag, New York, pages 339-342
*/

new;
library tsm,optmum,pgraph;

y = rndn(1000,1);

/*   Use the Kolmogorov-Smirnov test with a non white noise process
     For example:

     y = recserar(rndn(1000,1),0,0.5);
*/



_fourier = 0;
{lambda,I} = PDGM(y);
N = rows(I);
q = N/2;

I_ = I[1:q];
Yi = cumsumc(I_)/sumc(I_);
x = seqa(1,1,q);

k_alpha = -1.36~1.36;                          /* 5% level */
C_bound = (x-1)/(q-1)+k_alpha*(q-1)^(-0.5);

graphset;
  title("Kolmogorov-Smirnov test using the spectral analysis");
  _pdate = ""; _pnum = 2;
  _pltype = 6|1|1;
  scale(0|500,0|1);
  xy(x,Yi~C_bound);


