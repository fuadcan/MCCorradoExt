/*
** Procedure ClusOut:
**
** Purpose:
** Generates output for clustering algorithm.
**
** Format:
** ClusOut( C , ClusName );
**
** Input:
** C        Kx1 vector, containing number of cluster in which series are
**              included (returned from CLUSTER).
** ClusName Kx1 character vector, containing names of series.
**
** Output:
** None, the procedure prints output. If the output is to be captured in an
** ASCII then open an output-file before calling this procedure.
*/
proc(0) = ClusOut( C , ClusName );

local K , NoSer , NoClust , i , Omat, v, cp, num;

K = Rows(C);
NoClust = MaxC( C );
NoSer = SeqA(1,1,NoClust)~SumC( C.==(SeqA(1,1,NoClust)') );
NoSer = Rev(SortC(NoSer,2));

num=SumC( C.==(SeqA(1,1,NoClust)'));

"OUTPUT CLUSTER ALGORITHM:";
" ";
Format /ld 1,0;
"# Series     : " K ;
"# Clusters   : " NoClust ;
" ";
" ";
i = 0;
do while ( i < NoClust );
    i = i+1;

    Format /ld 1,0;
    "CLUSTER " i ;
    "Number of series in cluster : " NoSer[i,2];
    "Series in cluster: ";
    Omat = DelIf( ClusName , (C./=NoSer[i,1]) );
    Format /ld 8,0;
    $Omat';
    " ";
    " ";
endo;



/*print num;*/

v=(SeqA(1,1,MaxC(num)));
cp=counts(num,v);

print "Size Number";
print v~cp;



"END OF OUTPUT" ;
" ";
endp;
