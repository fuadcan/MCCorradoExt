new;
library tsm,optmum,pgraph;

rndseed 123456;

beta = -0.25;
Nobs = 500;

/* Fractional process with d = -0.25 */
{y1,retcode} = RND_arfima(beta,0,0,1,1000,Nobs,1);

/* Fractional process with d = 0.75 */
y2 = recserar(y1,0,1);

t = seqa(1,1,Nobs);

graphset;
  _pdate = ""; _pnum = 2; _pltype = 6;
  title("Fractional processes");
  _plegctl = {2 7 1 1}; _plegstr = "d=-0.25\0d=0.75";
  xlabel("time");
  xy(t,y1~y2);



