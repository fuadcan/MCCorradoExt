new;
library tsm,optmum,pgraph;
TSMset;

load purse[] = purse.asc;
Nobs = rows(purse);

/*
** State Space Model analysis
*/

{beta,stderr,Mcov,Logl} =  sm_LL(purse,1|1);

Z = 1; d_ = 0; H_ = beta[2]^2;
T = 1; c = 0; R = 1; Q = beta[1]^2;

call SSM_build(Z,d_,H_,T,c,R,Q,0);

a_0 = purse[1];
P_0 = zeros(1,1);

call KFiltering(purse,a_0,P_0);

a_SSM = KF_matrix(3);

/*
** GFLS analysis
*/

proc H(t);
  retp(1);
endp;

proc m(t);
  retp(0);
endp;

proc MM(t);
  retp(1);
endp;

proc F(t);
  retp(1);
endp;

proc d(t);
  retp(0);
endp;

proc DD(t);
  retp(1);
endp;

Q0 = 1;
p0 = purse[1];

{a1,u,as,us} = gfls(purse,&H,&m,&MM,&F,&d,&DD,Q0,p0,1);
{a10,u,as,us} = gfls(purse,&H,&m,&MM,&F,&d,&DD,Q0,p0,10);
{a100,u,as,us} = gfls(purse,&H,&m,&MM,&F,&d,&DD,Q0,p0,100);

graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _ptitlht = 0.30; _paxht = 0.40;
  fonts("simplex simgrma");
  begwind;
  window(1,2,0);
  setwind(1);
    title("LOCAL LEVEL MODEL"\
          "\Ly]t[ = \202b\201]t[+\202e\201]t[   "\
          "\202b\201]t[ = \202b\201]t-1[+\202c\201]t[");
    ylabel("\202b\201]t[");
    xy(seqa(1,1,Nobs),a_SSM);
  setwind(2);
    title("GFLS of the model"\
          "\Ly]t[ \202? b\201]t[   "\
          "\202b\201]t[ \202? b\201]t-1[");
    ylabel("\202b\201]t[");
    _plegstr = "\202m\201=1\0\202m\201=10\0\202m\201=100";
    _plegctl = {2 5 6 4.5};
    xy(seqa(1,1,Nobs),a1~a10~a100);
  endwind;

