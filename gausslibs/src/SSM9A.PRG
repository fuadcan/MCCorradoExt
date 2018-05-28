new;
library tsm,optmum,pgraph;
TSMset;

Z = {1 1,4 2,2 -3}; d = 0|0|0; H = {5 1 0,1 4 0,0 0 8};
T = {0.5 0.45,-0.5 0.8};
c = {0,0}; R = eye(2); Q = {2 0.5,0.5 1};

call SSM_build(Z,d,H,T,c,R,Q,0);

Nr = 20;
s = seqa(0,1,Nr+1);

e = 1|0;
{Delta,Ksi,zeta} = SSM_orthogonal(e,Nr);

graphset;
  begwind;
  makewind(9,0.5,0,6.855-0.5,0);
  makewind(9/2,5.855/2,0,6.855/2,0);
  makewind(9/2,5.855/2,4.5,6.855/2,0);
  makewind(9,0.5,0,5.855/2,0);
  makewind(9/2,5.855/2,0,0,0);
  makewind(9/2,5.855/2,4.5,0,0);

  setwind(1);
   title("Responses to the orthogonal impulse e = @[1 0@]'");
    _paxes = 0;
    _pnum = 0;
    _ptitlht = 2;
   draw;
  setwind(2);
    graphset;
    _ptitlht = 0.35; _paxht = 0.25; _pnumht = 0.25;
    _pnum = 2;
    _pline = 1~1~0~0~Nr~0~1~7~0;
    title("of the first variable y]1[(t)");
    xtics(0,Nr,1,0);
    xy(s,Delta[.,1]);
  setwind(3);
    title("of the second variable y]2[(t)");
    xtics(0,Nr,1,0);
    xy(s,Delta[.,2]);
  setwind(4);
   graphset;
   title("Cumulative responses");
    _paxes = 0;
    _pnum = 0;
    _ptitlht = 2;
   draw;
  setwind(5);
    graphset;
    _ptitlht = 0.35; _paxht = 0.25; _pnumht = 0.25;
    _pnum = 2;
    title("of the first variable y]1[(t)");
    xtics(0,Nr,1,0);
    _pline = 1~1~0~zeta[1]~Nr~zeta[1]~1~7~0;
    xy(s,Ksi[.,1]);
  setwind(6);
    title("of the second variable y]2[(t)");
    xtics(0,Nr,1,0);
    _pline = 1~1~0~zeta[2]~Nr~zeta[2]~1~7~0;
    xy(s,Ksi[.,2]);
  endwind;





