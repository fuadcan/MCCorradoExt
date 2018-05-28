/*
** Procedure MVKPSS:
**
** Purpose:
** Calculates the multivariate KPSS-test.
**
** Format:
** { W0 , W1 } = MVKPSS( Y , l );
**
** Input:
** Y    (TxN)-matrix, containing series.
** l    scalar, truncation parameter for Bartlett window.
**
** Output:
** W0   scalar, value of multivariate KPSS test for null hypothesis of zero
**          mean stationarity.
** W1   scalar, value of multivariate KPSS test for null hypothesis of level
**          stationarity.
*/
Proc(2) = MVKPSS( Y , l );



local psi1,T,S,Sdemean,W0,W1,i,X,B,e;

T = rows(Y);

X = ones(T,1)~seqa(1,1,T);
B = Y/X;
e = Y-X*B;

psi1 = NeweyWst( e , l );


psi1 = invpd(psi1);

S = cumsumc(Y);
Sdemean = cumsumc( (Y-meanc(Y)') );

W0 = 0;
W1 = 0;

i = 0;
do while (i < T);
    i = i + 1;

        locate 9,1;
        format /ld 8,0;
        
        

    W0 = W0+S[i,.]*psi1*S[i,.]';
    W1 = W1+Sdemean[i,.]*psi1*Sdemean[i,.]';


endo;

W0 = W0./(T^2);
W1 = W1./(T^2);

retp( W0 , W1 );



endp;


