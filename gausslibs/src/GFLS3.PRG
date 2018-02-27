new;
library pgraph,tsm;
TSMset;

rndseed 123;

Nobs = 500;
t = seqa(1,1,Nobs);
x = floor(rndu(Nobs,3)*100);

S = 250;

b1 = ones(S,1)|ones(Nobs-S,1)*2;
b2 = recserar(rndn(Nobs,1),0,.5);
b3 = seqa(1,10/Nobs,Nobs);
B = b1~b2~b3;

y1 = x[.,1].*b1 + x[.,3].*b3 + rndn(Nobs,1);
y2 = x[.,2].*b2 + x[.,3].*b3 + rndn(Nobs,1);

proc H(t);
  local Ht;
  Ht = (x[t,1]~0~x[t,3])|
       (0~x[t,2]~x[t,3]);
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
p0 = zeros(3,1);

{a,u,as,us} = gfls(y1~y2,&H,&m,&MM,&F,&d,&DD,Q0,p0,1);

proc DD2(t);
  local Dt;
  Dt = eye(3)/10; Dt[1,1] = 1;
  retp(Dt);
endp;

{a2,u2,as2,us2} = gfls(y1~y2,&H,&m,&MM,&F,&d,&DD2,Q0,p0,10);

i = 1;

do until i > 3;

  graphset;
    _pdate = "";

    begwind;
    makewind(9,1,0,5.855,0);
    makewind(4.5,5.855/2,0,5.855/2,0);
    makewind(4.5,5.855/2,4.5,5.855/2,0);
    makewind(4.5,5.855/2,0,0,0);
    makewind(4.5,5.855/2,4.5,0,0);

    setwind(1);
      _paxes = 0;
      _pnum = 0;
      _ptitlht = 1;

      str = ftos(i,"Parameter #%lf",1,0);
      title(str);
      draw;

    graphset;
    _pdate = ""; _pnum = 2; _pnumht = 0.25; _paxht = 0.30; _ptitlht = 0.30;

    setwind(2);
      title("True values");
      xy(t,B[.,i]);
    setwind(3);
      title("GFLS estimate (first model)");
      xy(t,a[.,i]);
    setwind(4);
      title("GFLS estimate (second model)");
      xy(t,a2[.,i]);
    setwind(5);
      title("Difference between the true \Lvalues and the GFLS estimate");
      xy(t,B[.,i]-(a[.,i]~a2[.,i]));
    endwind;

  i = i + 1;

endo;
