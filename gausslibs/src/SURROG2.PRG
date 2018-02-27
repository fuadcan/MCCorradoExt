new;
library tsm,pgraph,optmum;

Nobs = 200;
s = seqa(1,1,Nobs);

let SIGMA = {1 0.5,0.5 1};
Y0 = 25|-12;
x = RND_arma(0.6|0.3|0|0.5,1,0,SIGMA,Y0,Nobs);
y = surrogate(x);

graphset;
  _pdate = ""; _pnum = 2;
  title("Surrogate data technique\Lfirst series");
  xlabel("Observation");
  _plegstr = "Original data\0Surrogate data";
  _plegctl = 1;
  xy(s,x[.,1]~y[.,1]);

graphset;
  _pdate = ""; _pnum = 2;
  title("Surrogate data technique\Lsecond series");
  xlabel("Observation");
  _plegstr = "Original data\0Surrogate data";
  _plegctl = 1;
  xy(s,x[.,2]~y[.,2]);

i = 0;
do until i>2;
  COVx = autoc(x,i);
  COVy = autoc(y,i);
  print ftos(i,"Autocovariance matrix of order %lf",1,0);
  print "          original data                 surrogate data";
  print COVx~COVy; print ;
  i = i +1;
endo;



proc autoc(x,i);
  local K,xi,data,retcode,cov;
  K = cols(x);
  xi = lagn(x,i);
  {data,retcode} = Missing(x~xi,0);
  x = data[.,1:K];
  xi = data[.,1+K:2*K];
  cov = x'xi/rows(x);
  retp(cov);
endp;

