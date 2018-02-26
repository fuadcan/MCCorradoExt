/*
** Procedure PValue:
**
** Purpose:
** Calculates the p-value of the multivariate KPSS test statistic.
**
** Format:
** { p0 , p1 } = pvalue( W0 , W1 , k , asymdis0 , asymdis1 );
**
** Input:
** W0       scalar, value of multivariate KPSS test for null hypothesis of zero
**              mean stationarity.
** W1       scalar, value of multivariate KPSS test for null hypothesis of level
**              stationarity.
** k        scalar, dimension of vector for which test values are determined.
** asymdis0 (lx99)-matrix, containing the first 99 percentiles of the asymptotic
**              distribution of test for zero mean stationarity in each column.
**              The row number corresponds to the dimension of the test.
** asymdis1 (lx99)-matrix, containing the first 99 percentiles of the asymptotic
**              distribution of test for zero mean stationarity in each column.
**              The row number corresponds to the dimension of the test.
**
** Output:
** p0       scalar, p-value associated with w0.
** p1       scalar, p-value associated with w1.
**
*/
proc(2) = pvalue( w0 , w1 , k , W0PH , W1PH );

local l , yint , xint0 , xint1 , p0 , p1 , m0 , m1 , ind0 , ind1 ;

_poldeg = 10 ;

l = rows(W0PH);
m0 = cols(W1PH);
m1 = cols(W1PH);

/* if ( k < 0 ) or ( k > l );
    Errorlog("Cannot determine p-value for test of this dimension");
    end;
endif; */

/* xint0 = asymdis0[k,.]';
xint1 = asymdis1[k,.]'; */

xint0 = W0PH;
xint1 = W1PH;

/*CHECKS NUMERIC VECTOR AGAINST ANOTHER and RETURNS THE INDICES OF ELEMENTS OF FIRST 
 VECTOR  IN SECOND VECTOR. IN OUR CASE AS W0PH and W1PH 
ARE SCALARS THIS CHECKS WHETHER W0PH>WO. if THIS IS TRUE IT RETURNS THE
INDEX OF THE ELEMENTS for WHICH THE CONDITION IS SATISFIED. AS HERE WE HAVE JUST ONE
SCALAR IT RETURNS 1 THAT IS INDEX OF TE ONLY ELEMENT WE HAVE GOT*/

ind0  = indnv( 1 , (xint0.>w0) ); /*CHECKS WHETHER WO<W0PH I.E. if WE ARE INSIDE THE
                                                        CRITICAL REGION or NOT*/
ind1  = indnv( 1 , (xint1.>w1) );


if ( ind0 == 1 );
    p0 = 0;
elseif ismiss( ind0 ); 
/*THIS RETURNS ONE if A MATRIX HAS A MISSING VALUE. IN FACT 
INDNV MAY return MISSING VALUE if THE CONDITION W0PH>WO IS NOT SATISFIED*/

    p0 = 1;
else;
    p0 = (m0+1-ind0)./(m0+1); /*THIS PART ONLY APPLIES if YOU HAVE ALL THE 
PERCENTILES OF THE ASYMPTOTYC DISTRIBUTION. IN OUR CASE WE HAVE JUST ONE SCALAR*/

    p0 = p0+(1./(m0+1)).*((w0-xint0[m0-1])./(xint0[m0]-xint0[m0-1]));
endif;

if ( ind1 == 1 );
    p1 = 0;
elseif ismiss( ind1 );
    p1 = 1;
else;
    p1 = (m1+1-ind1)./(m1+1);
    p1 = p1+(1./(m1+1)).*((w1-xint1[m1-1])./(xint1[m1]-xint1[m1-1]));
endif;

p0 = 1 - p0;
p1 = 1 - p1;

retp( p0 , p1 );
endp;


