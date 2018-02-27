new;
library tsm,optmum,pgraph;
load gnp[] = gnp.asc;
Nobs = rows(gnp);


/*
** State Space Model analysis
*/

{beta,stderr,Mcov,Logl} =  sm_LLT(gnp,1|1|1);

Z = 1~0; d_ = 0;
T = {1 1,0 1}; c = 0|0; R = {1 0,0 1};
theta = beta^2;  /* sigma^2 */
H_ = theta[3];
Q = (theta[1]~0)|(0~theta[2]);

call SSM_build(Z,d_,H_,T,c,R,Q,0);

a0 = gnp[1]|0; P0 = zeros(2,2);
call KFiltering(gnp,a0,P0);

a_SSM = KF_matrix(3);

/*
** GFLS analysis
*/

proc H(t);
  local Ht;
  Ht = (1~0);
  retp(Ht);
endp;

proc m(t);
  retp(0);
endp;

proc MM(t);
  retp(1);
endp;

proc F(t);
  local Ft;
  Ft = (1~1)|(0~1);
  retp(Ft);
endp;

proc d(t);
  retp(0|0);
endp;

proc DD(t);
  retp(eye(2));
endp;

Q0 = eye(2);
p0 = gnp[1]|0;

{a1,u,as,us} = gfls(gnp,&H,&m,&MM,&F,&d,&DD,Q0,p0,1);
{a10,u,as,us} = gfls(gnp,&H,&m,&MM,&F,&d,&DD,Q0,p0,10);
{a100,u,as,us} = gfls(gnp,&H,&m,&MM,&F,&d,&DD,Q0,p0,100);

graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _ptitlht = 0.30; _paxht = 0.40;
  fonts("simplex simgrma");
  begwind;
  window(2,1,0);
  setwind(1);
    title("LOCAL LINEAR TREND MODEL"\
          "\Ly]t[=\202d\201]t[+\202e\201]t[  "\
          "\202d\201]t[=\202d\201]t-1[+\202b\201]t-1[+\202c\201]t[  "\
          "\202b\201]t[=\202b\201]t-1[+\202z\201]t[");
    ylabel("\202d\201]t[");
    xy(seqa(1,1,Nobs),a_SSM[.,1]);
  setwind(2);
    title("GFLS of the model"\
          "\Ly]t[ \202? d\201]t[  "\
          "\202d\201]t[ \202? d\201]t-1[+\202b\201]t-1[  "\
          "\202b\201]t[ \202? b\201]t-1[");
    ylabel("\202d\201]t[");
    _plegstr = "\202m\201=1\0\202m\201=10\0\202m\201=100";
    _plegctl = {2 5 7 1};
    xy(seqa(1,1,Nobs),a1[.,1]~a10[.,1]~a100[.,1]);
  endwind;


graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.25; _ptitlht = 0.30; _paxht = 0.40;
  fonts("simplex simgrma");
  begwind;
  window(2,1,0);
  setwind(1);
    title("LOCAL LINEAR TREND MODEL"\
          "\Ly]t[=\202d\201]t[+\202e\201]t[  "\
          "\202d\201]t[=\202d\201]t-1[+\202b\201]t-1[+\202c\201]t[  "\
          "\202b\201]t[=\202b\201]t-1[+\202z\201]t[");
    ylabel("\202b\201]t[");
    xy(seqa(1,1,Nobs),a_SSM[.,2]);
  setwind(2);
    title("GFLS of the model"\
          "\Ly]t[ \202? d\201]t[  "\
          "\202d\201]t[ \202? d\201]t-1[+\202b\201]t-1[  "\
          "\202b\201]t[ \202? b\201]t-1[");
    ylabel("\202b\201]t[");
    _plegstr = "\202m\201=1\0\202m\201=10\0\202m\201=100";
    _plegctl = {2 5 7 1};
    xy(seqa(1,1,Nobs),a1[.,2]~a10[.,2]~a100[.,2]);
  endwind;



