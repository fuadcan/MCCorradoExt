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
proc(2) = pvalue( W0 , W1 , k , W0PH , W1PH );
local l , yint , xint0 , xint1 , p0 , p1 , m0 , m1 , ind0 , ind1;

_poldeg = 10 ;

l = cols(W0PH);

m0 = rows(W0PH);

m1 = rows(W1PH);





/*l = rows(W0PH);   we are commenting this out 
m0 = cols(W1PH);
m1 = cols(W1PH); */

/* if ( k < 0 ) or ( k > l );
    Errorlog("Cannot determine p-value for test of this dimension");
    end;
endif; */

/* xint0 = asymdis0[k,.]';
xint1 = asymdis1[k,.]'; */

xint0 = W0PH;
xint1 = W1PH;

/*
- THE FUNCTION INDNV NUMERIC VECTOR AGAINST ANOTHER and RETURNS THE INDICES OF 
ELEMENTS OF FIRST  VECTOR (OR SCALAR)  IN SECOND VECTOR. 

 - IN OUR CASE THE FUNCTION RETURNS THE
INDEX OF THE FIRST ELEMENT FOR WHICH THE CONDITION (xint0.>w0) IS SATISFIED
(I.E. WHERE IT FIRST FINDS THE SCALAR 1 IN THE VECTOR (xint0.>w0)*/


ind0  = indnv( 1 , (xint0.>W0) ); 

/*
- it checks for each element whether  w0 exceeds critical value i.e. if we are 
inside the critical region or not.

- the function indnv returns the index of the (first) percentile for the k-th row
where condition is true */
ind1  = indnv( 1 , (xint1.>W1) );



if ( ind0 == 1 ); /* if the index 
                  of the percentile for which the condition (xint0.>w0) is 1
                  then it sets p0=0*/

    p0 = 0;
elseif ismiss( ind0 ); /*THIS RETURNS p0=1 if ind0 HAS A MISSING VALUE.
 IN FACT INDNV MAY return MISSING VALUE if THE CONDITION W0PH>WO IS NOT SATISFIED*/

    p0 = 1;
else;
    p0 = (m0+1-ind0)./(m0+1);  /*THIS PART ONLY APPLIES if YOU HAVE that the index 
                               of THE  PERCENTILES for which for which 
                           the condition (xint0.>w0) is different from zero or 1*/

    p0 = p0+(1./(m0+1)).*((W0-xint0[m0-1])./(xint0[m0]-xint0[m0-1]));


endif;


if ( ind1 == 1 );
    p1 = 0;
elseif ismiss( ind1 );
    p1 = 1;
else;
    p1 = (m1+1-ind1)./(m1+1);
    p1 = p1+(1./(m1+1)).*((W1-xint1[m1-1])./(xint1[m1]-xint1[m1-1]));

endif;

p0 = 1 - p0;
p1 = 1 - p1;

retp( p0 , p1 );

endp;


