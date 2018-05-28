=======================================================================
ASYMPTOTICALLY PERFECT AND RELATIVE CONVERGENCE OF PRODUCTIVITY
Programs, August 1998. Bart Hobijn and Philip Hans Franses
=======================================================================
The programs and data consist of the following files:

The following files are in CONVERGE.ZIP:

Example programs:
PWT      GSS         1,812  12-15-98  6:32p PWT.GSS
BERNDURL GSS         1,859  12-15-98  6:29p BERNDURL.GSS

GAUSS-procedures:
LOCMIN   G             447  08-02-98 10:11p LOCMIN.G
COINTMAT G             479  08-29-98 12:01a COINTMAT.G
CLUSTER  G           7,704  12-15-98  6:18p cluster.g
CLUSOUT  G           1,145  03-12-95  1:25p CLUSOUT.G
CLUSCORR G             669  08-20-98  2:38p CLUSCORR.G
MVKPSS   G             951  08-21-98 12:23p MVKPSS.G
PVALUE   G           2,000  08-20-98 10:33a PVALUE.G
NEWEYWST G             782  07-24-98  3:00p NEWEYWST.G
BARTLETT G             599  01-01-97  7:51p BARTLETT.G
LAGMAT   G             103  07-22-98  9:29p LAGMAT.G

Datafiles to calculate p-values (used by cluster.g):
ASYMDST0 DAT        39,600  08-03-98  5:49a ASYMDST0.DAT
ASYMDST0 DHT           920  08-03-98  5:49a ASYMDST0.DHT
ASYMDST1 DAT        39,600  08-03-98  5:49a ASYMDST1.DAT
ASYMDST1 DHT           920  08-03-98  5:49a ASYMDST1.DHT

Datafiles:
PWT      TXT        16,908  04-16-99 11:52a PWT.TXT    (ASCII)
PWT      DAT        26,880  12-15-98  6:31p PWT.DAT    (GAUSS)
PWT      DHT         1,024  12-15-98  6:31p PWT.DHT    (GAUSS header file)
BERNDURL TXT         7,188  04-16-99 11:54a BERNDURL.TXT (ASCII)
BERNDURL DHT           256  03-19-97  7:19p BERNDURL.DHT (GAUSS)
BERNDURL DAT        11,264  03-19-97  7:19p BERNDURL.DAT (GAUSS header file)

This file:
README   TXT                12-15-98  6:37p readme.txt

============================================================================
Datafiles:

Both data files contain per capita RGDP levels in ppp dollars with baseyear
1980 for the Bernard and Durlauf-data and 1985 for the PWT-data.

The datafiles have the following format.
=>  Series are in columns.
=>  Headers contain the country abbreviation corresponding to the series in
    the respective column.
=>  First column contains the year.

Sources:
PWT data                    taken from the Penn World Table, Mark 5.5
Bernard and Durlauf data    taken from Bernard, A.B. and S.N. Durlauf (1995),
                            "Convergence in International Output",
                            Journal of Applied Econometrics, 10, 161-173.
============================================================================
GAUSS-procedures:
Hopefully the comments suffice when using them. They have been tested in
GAUSS 3.2.1 for DOS. If you have any questions or problems, please contact
Bart Hobijn (hobijnb@fasecon.econ.nyu.edu).
The main procedure is CLUSTER.G. The example programs illustrate how to use it.
