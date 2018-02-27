new;
library tsm,optmum;
TSMset;

rndseed 123456;

et = rndn(100,1)*1;
et_1 = shiftr(et',1,0)';
bt = recserar(et-0.5*et_1,10,1);      /* ARMA(1,1) process */

xt = rndu(100,1);
yt = xt.*bt + rndn(100,1)*2;

save kal3a = yt;
save kal3b = xt;

