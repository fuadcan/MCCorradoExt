new;
library pgraph,tsm;
TSMset;

_print = 0;

rndseed 123;

Nobs = 1000;
t = seqa(1,1,Nobs);

sy = 8*cos(seqa(0,2*pi/Nobs,Nobs)) + rndn(Nobs,1)*0.25;
sx = 10*sin(seqa(0,2*pi/Nobs,Nobs)) + rndn(Nobs,1)*0.25;
c = 2*recserar(rndn(Nobs,1),10,0.9);

y = sy + c + rndn(Nobs,1);
x = sx + c + rndn(Nobs,1);

proc H(t);
  local Ht;
  Ht = (1~0~1)|
       (0~1~1);
  retp(Ht);
endp;

proc m(t);
  retp(0|0);
endp;

proc MM(t);
  retp(eye(2));
endp;

proc F(t);
  retp(eye(3));
endp;

proc d(t);
  retp(zeros(3,1));
endp;

proc DD(t);
  retp(eye(3));
endp;

Q0 = eye(3);
p0 = 0|0|0;
{a,u,as,us} = gfls(y~x,&H,&m,&MM,&F,&d,&DD,Q0,p0,1);

p0 = y[1]|x[1]|0;
{a2,u2,as2,us2} = gfls(y~x,&H,&m,&MM,&F,&d,&DD,Q0,p0,1);

graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.15; _ptitlht = 0.15; _ptitlht = 0.15;
  title("specific component of the y series");
  _pltype = 6|1|3;
  _plegstr ="true\0GFLS (first model)\0GFLS (second model)";
  _plegctl = {2 4 3 5};
  xy(t,sy~as[.,1]~as2[.,1]);

  title("specific component of the x series");
  _plegctl = {2 4 2 0.75};
  xy(t,sx~as[.,2]~as2[.,2]);

  title("common component");
  _plegctl = {2 4 1 0.75};
  xy(t,c~as[.,3]~as2[.,3]);



