/*
**  Prichard, D. and J. Theiler [1994], Generating surrogate data for
**  time series with several simultaneously measured variables,
**  working paper, Santa Fe Institute
**  (http://www.santafe.edu/...)
**
**  Palus, M. [1995], Detecting nonlinearity in multivariate time series,
**  working paper, Santa Fe Institute
**  (http://www.santafe.edu/...)
**
**  Lorenz model
**
*/

new;
library tsm,optmum,pgraph;
#include rk4.src;
TSMset;

declare external sigma,beta,R;

proc Lorenz(t,u);
  local x,y,z,dx,dy,dz,du;

  x = u[1]; y = u[2]; z = u[3];


  dx = sigma*(y-x);
  dy = -x*z + R*x - y;
  dz = x*y - beta*z;

  du = dx|dy|dz;

  retp(du);
endp;


/*
**  Initializing of the parameters and the starting
**  position of the system
*/

sigma = 10; beta = 8/3; R = 28;
u0 = 0|1|1;

/*
{t,u} = rk4(&Lorenz,u0,0,5.5,512);
u0 = u[512,.]';
*/
u0 = 15.34|13.68|37.91;
{t,u} = rk4(&Lorenz,u0,0,40,1028);

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

    title("Lorenz Model");
    draw;

  graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _paxht = 0.30;

  setwind(2);
    xlabel("x]t[");
    ylabel("y]t[");
    zlabel("z]t[");
    xyz(u[.,1],u[.,2],u[.,3]);

  setwind(3);
    xlabel("x]t[");
    ylabel("y]t[");
    xy(u[.,1],u[.,2]);

  setwind(4);
    xlabel("x]t[");
    ylabel("z]t[");
    xy(u[.,1],u[.,3]);

  setwind(5);
    xlabel("y]t[");
    ylabel("z]t[");
    xy(u[.,2],u[.,3]);

  endwind;


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

    title("Lorenz Model");
    draw;

  graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _paxht = 0.30;

  setwind(2);
    xlabel("x]t[");
    ylabel("x]t-1[");
    xy(u[.,1],lagn(u[.,1],1));

  setwind(3);
    xlabel("x]t[");
    ylabel("x]t-2[");
    xy(u[.,1],lagn(u[.,1],2));

  setwind(4);
    xlabel("x]t[");
    ylabel("x]t-5[");
    xy(u[.,1],lagn(u[.,1],5));

  setwind(5);
    xlabel("x]t[");
    ylabel("x]t-10[");
    xy(u[.,1],lagn(u[.,1],10));

  endwind;




_fourier = 0;
v = surrogate(u);


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

    title("FT surrogate data of the Lorenz Model");
    draw;

  graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _paxht = 0.30;

  setwind(2);
    xlabel("x]t[");
    ylabel("y]t[");
    zlabel("z]t[");
    xyz(v[.,1],v[.,2],v[.,3]);

  setwind(3);
    xlabel("x]t[");
    ylabel("y]t[");
    xy(v[.,1],v[.,2]);

  setwind(4);
    xlabel("x]t[");
    ylabel("z]t[");
    xy(v[.,1],v[.,3]);

  setwind(5);
    xlabel("y]t[");
    ylabel("z]t[");
    xy(v[.,2],v[.,3]);

  endwind;


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

    title("FT surrogate data of the Lorenz Model");
    draw;

  graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _paxht = 0.30;

  setwind(2);
    xlabel("x]t[");
    ylabel("x]t-1[");
    xy(v[.,1],lagn(v[.,1],1));

  setwind(3);
    xlabel("x]t[");
    ylabel("x]t-2[");
    xy(v[.,1],lagn(v[.,1],2));

  setwind(4);
    xlabel("x]t[");
    ylabel("x]t-5[");
    xy(v[.,1],lagn(v[.,1],5));

  setwind(5);
    xlabel("x]t[");
    ylabel("x]t-10[");
    xy(v[.,1],lagn(v[.,1],10));

  endwind;


