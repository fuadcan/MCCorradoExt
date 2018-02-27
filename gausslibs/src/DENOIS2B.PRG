/*
** Donoho and Johnstone [1994], Ideal spatial adaptation by wavelet shrinkage,
** Biometrika, 81, 425-455
*/

new;
library tsm,pgraph;

rndseed 123456;

proc sgn(x);
  local s;
  s = 2*(x.>0)-1;
  retp(s);
endp;

proc K1(t);
  retp(0.5+0.5*sgn(t));
endp;

proc K2(t);
  retp((1+abs(t))^(-4));
endp;

epsilon = 0.05;

N = 2048;
t = seqa(0,1/(N-1),N);

t1 = 0.1|0.13|0.15|0.23|0.25|0.40|0.44|0.65|0.76|0.78|0.81;
h1 = 4|-5|3|-4|5|-4.2|2.1|4.3|-3.1|2.1|-4.2;
t2 = t1;
h2 = 4|5|3|4|5|4.2|2.1|4.3|3.1|5.1|4.2;
w2 = 5|6|6|10|10|30|10|10|5|8|5; w2 = w2/1000;


Blocks = 4*sumc(h1.*K1(t'-t1));
Bumps = 10*sumc(h2.*K2((t'-t2)./w2));
HeaviSine = 2.5*4*sin(4*pi*t) - sgn(t-0.3) - sgn(0.72-t);
Doppler = 10*sqrt(t.*(1-t)).*sin(2*pi*(1+epsilon)./(t+epsilon));

graphset;
begwind;
window(2,2,0);
_pdate = ""; _ptitlht = 0.30; _pnumht = 0.20;
setwind(1);
  title("Blocks");
  xy(t,Blocks);
nextwind;
  title("Bumps");
  xy(t,Bumps);
nextwind;
  title("HeaviSine");
  xy(t,HeaviSine);
nextwind;
  title("Doppler");
 xy(t,Doppler);
endwind;


sc = stdc(Blocks)/14;
Noisy_Blocks = Blocks + rndn(N,1)*sc;
sc = stdc(Bumps)/14;
Noisy_Bumps = Bumps + rndn(N,1)*sc;
sc = stdc(HeaviSine)/14;
Noisy_HeaviSine = HeaviSine + rndn(N,1)*sc;
sc = stdc(Doppler)/14;
Noisy_Doppler = Doppler + rndn(N,1)*sc;

graphset;
begwind;
window(2,2,0);
_pdate = ""; _ptitlht = 0.30; _pnumht = 0.20;
setwind(1);
  title("Noisy Blocks");
  xy(t,Noisy_Blocks);
nextwind;
  title("Noisy Bumps");
  xy(t,Noisy_Bumps);
nextwind;
  title("Noisy HeaviSine");
  xy(t,Noisy_HeaviSine);
nextwind;
  title("Noisy Doppler");
 xy(t,Noisy_Doppler);
endwind;


{H,G,Htilde,Gtilde} = Daubechies(12);

_wcenter = 0;
x = Noisy_Blocks~Noisy_Bumps~Noisy_HeaviSine~Noisy_Doppler;
y = zeros(N,4);

i = 1;
do until i>4;

  pkt = wpkt(x[.,i],H,G,0);
  B = 1|2|3|4|4;
  wpk = Basis(pkt,B);
  bnd1 = wpk[1:N/2];
  bnd2 = wpk[N/2+1:N];
  bnd2 = Thresholding(bnd2,0.99);
  wpk = bnd1|bnd2;

  y[.,i] = iwpkt(wpk,Htilde,Gtilde,B);

  i = i + 1;
endo;


graphset;
begwind;
window(2,2,0);
_pdate = ""; _ptitlht = 0.30; _pnumht = 0.20;
setwind(1);
  title("Blocks");
  xy(t,y[.,1]);
nextwind;
  title("Bumps");
  xy(t,y[.,2]);
nextwind;
  title("HeaviSine");
  xy(t,y[.,3]);
nextwind;
  title("Doppler");
 xy(t,y[.,4]);
endwind;



