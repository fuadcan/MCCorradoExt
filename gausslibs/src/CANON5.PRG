/*
**  Compute the Impulse Responses and the Accumulated Impulse Responses
**  of an ARMA process
*/

new;
library tsm,optmum,pgraph;

Nresp = 10;
Period = seqa(1,1,Nresp);

beta = 0.5|-0.15|0.4|0.4|-0.3|0;
{phi,theta,retcode} = canonical_arfima(beta,3,2,0,Nresp-1);

graphset;
  _pdate = ""; _pnum = 2; _pnumht = 0.15; _ptitlht = 0.15;
    title("Responses of the innovations of the ARMA process"\
          "\Ly]t[-0.5y]t-1[+0.15y]t-2[-0.4y]t-3[ = "\
          "u]t[-0.4u]t-1[+0.3u]t-2[");
    xlabel("periods");
    xtics(1,Nresp,1,0);
    _plegstr = "Impulse Responses\0Accumulated Responses";
    _plegctl = {2 5 4 2};
    xy(Period,theta~cumsumc(theta));
