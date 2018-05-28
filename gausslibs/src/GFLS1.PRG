new;
library pgraph,tsm;

rndseed 123;
Nobs = 200;

x = floor(rndu(Nobs,2)*100);
x = ones(Nobs,1)~x;

b1 = ones(Nobs,1);
b2 = reshape(8|0|-3|-5,Nobs,1);
b2 = b2/maxc(abs(b2));
b3 = recserar(rndn(Nobs,1),10,0.9);
b3 = b3/maxc(abs(b3));
B = b1~b2~b3;

y = sumc(x'.*B') + rndn(Nobs,1);

B = eye(3);
C = eye(3);
D = eye(3).*(1000|0|1);
F = eye(3).*(0|1|0);

_print = 0;
{beta,u,SSE} = GFLS2(y,x,4,1|100,B,C,D,F);
{Bfls,u,r2M,r2D} = FLS(y,x,1000|0.01|1);


t = seqa(1,1,Nobs);

graphset;
  _pdate = ""; _pnum = 2;
  title("Generated time series");
  xy(t,y);

graphset;
  _pdate = ""; _pnum = 2;
  scale(0|nobs,-1|5); _paxht = 0.25; _ptitlht = 0.25; _pnumht = 0.25;
  begwind;
  window(1,2,0);
    title("First coefficient");
    xy(t,b1);
  nextwind;
    title("Estimated");
    _plegstr = "GFLS\0FLS";
    _plegctl = {2 7 2 5};
    xy(t,beta[.,1]~Bfls[.,1]);
  endwind;

graphset;
  _pdate = ""; _pnum = 2;
  scale(0|nobs,-1|2); _paxht = 0.25; _ptitlht = 0.25; _pnumht = 0.25;
  begwind;
  window(1,2,0);
    title("Second coefficient");
    xy(t,b2);
  nextwind;
    title("Estimated");
    _plegstr = "GFLS\0FLS";
    _plegctl = {2 7 2 5};
    xy(t,beta[.,2]~Bfls[.,2]);
  endwind;

graphset;
  _pdate = ""; _pnum = 2;
  scale(0|nobs,-1|2); _paxht = 0.25; _ptitlht = 0.25; _pnumht = 0.25;
  begwind;
  window(1,2,0);
    title("Third coefficient");
    xy(t,b3);
  nextwind;
    title("Estimated");
    _plegstr = "GFLS\0FLS";
    _plegctl = {2 7 2 5};
    xy(t,beta[.,3]~Bfls[.,3]);
  endwind;

