/*
** Procedure LocMin:
** Gives row and column of minimum element in matrix
** If minimum is not unique it returns the minimum element with the lowest
** column and if this column contains 2 elements equal to the minimum it
** returns the lowest row.
**
*/
Proc(2) = LocMin( X );

Local MinX , r , c ;

MinX = MinC(MinC(X));

c = IndNV( 1 , ( SumC( X.==MinX )./=0 ) );
r = IndNV( 1 , ( X[.,c].==MinX ) );

RetP( r , c );
EndP;
