new;
library tsm,pgraph,optmum;

Nobs = 500;
s = seqa(1,1,Nobs);

let SIGMA = {1 0.5,0.5 1};
Y0 = miss(0,0);
x = RND_arma(0.6|0.3|0|0.5,1,0,SIGMA,Y0,Nobs);
x = x - meanc(x)';

x1 = x[.,1]; x2 = x[.,2];

_fourier = 0;

{lambda,I} = PDGM(x1);
cov = real(inverse_fourier(I));

R1 = autoc(x1,10);

graphset;
    _pdate = "";
    title("Covariogram");
    _pltype = 6|1;
    xlabel("Lag");
    xy(seqa(0,1,11),R1~cov[1:11]);


proc autoc(x,k);
  local t,rho;
  x=packr(x); t=rows(x); x=x-meanc(x);
  rho=rev(conv(x,rev(x),t-k,t));
  retp(rho/t);
endp;

