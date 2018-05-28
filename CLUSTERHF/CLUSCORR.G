/*
** Procedure CLUSCORR:
**
** Purpose:
** Computes cluster-correlation.
**
** Format:
** r = ClusCorr( C1 , C2 );
**
** Input:
** C1       Nx1 vector, outcome of cluster procedure
** C2       Nx1 vector, outcome of cluster procedure
**
** Output:
** r        scalar, cluster-correlation between C1 and C2.
*/
Proc(1) = ClusCorr( C1 , C2 );

Local D1 , D2 , n , r;

n = Rows(C1);
If ( Rows(C2) /= n );
    ErrorLog("C1 and C2 must have same number of rows");
    End;
EndIf;

D1 = ( C1 .== C1' );
D2 = ( C2 .== C2' );

r = (SumC(SumC(D1.*D2))-n);
r = r./SqRt( (SumC(SumC(D1))-n).*(SumC(SumC(D2))-n) );
r = sqrt(r);

RetP(r);
EndP;

