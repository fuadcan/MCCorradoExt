/*
** Procedure NeweyWst:
**
** Purpose:
** Estimates sample periodogram at frequency zero, using the estimator
** introduced by Newey and West (1986,EMA).
**
** Format:
** S = NeweyWst( Y , l );
**
** Input:
** Y    (TxN)-matrix, containing series.
** l    scalar, truncation parameter for the Bartlett window.
**
** Output:
** S    (NxN)-matrix, containing estimated sample periodogram at freq. zero.
*/
Proc(1) = NeweyWst( Y , l );

local gami , T , S , i , Yl , kernel;

T = rows(Y);
Y = (Y-meanc(Y)');


S = (Y'*Y)./T;
Yl = Y;
kernel = Bartlett( T , l );

i = 0;
do while ( i < l );
    i = i+1;

    Yl = lagmat(Yl);
    gami = (Y[i+1:rows(Y),.]'*Yl[i+1:rows(Yl),.])./T;
    S = S+kernel[i+1].*(gami+gami');
endo;

retp(S);
endp;




