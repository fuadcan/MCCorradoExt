/*
** Procedure Bartlett:
**
** Purpose:
** Calculates Bartlett kernel.
**
** Format:
** kern = Bartlett( T , m )
**
** Input:
** T    scalar, sample size.
** m    scalar, bandwidth.
**
** Output:
** kern (Tx1)-vector, kernel weights for the estimated autocovariances for the
**       estimation of the smoothed periodogram.
**
** Note:
** k(m,j) = ( 1 - (j/(m+1)) )       for j <= m
**        = 0                       otherwise
*/
proc(1) = bartlett( T , m );

local j , kern;

j = seqa(0,1,T);
kern = ( 1 - j./(m+1));
kern = kern.*(kern.>0);

retp(kern);
endp;


