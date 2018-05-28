Bart Hobijn and Philip Hans Franses, "Asymptotically Perfect and Relative
Convergence of Productivity," Journal of Applied Econometrics, Vol. 15,
No. 1, 2000, pp. 59-81.

The programs and data consist of the following files, all of which are in
DOS format and zipped in the file hobijn-franses.zip:

Example programs:
PWT      GSS         1,812  12-15-98  6:32p pwt.gss
BERNDURL GSS         1,859  12-15-98  6:29p berndurl.gss
EURO2009 GSS         CHANGE DATAFILE FOR EUROPEAN SECTORAL GVA  

Gauss Procedures:
LOCMIN   G             447  08-02-98 10:11p locmin.g
COINTMAT G             479  08-29-98 12:01a cointmat.g
CLUSTER  G           7,704  12-15-98  6:18p cluster.g
CLUSOUT  G           1,145  03-12-95  1:25p clusout.g
CLUSCORR G             669  08-20-98  2:38p cluscorr.g
MVKPSS   G             951  08-21-98 12:23p mvkpss.g
PVALUE   G           2,000  08-20-98 10:33a pvalue.g
NEWEYWST G             782  07-24-98  3:00p neweywst.g
BARTLETT G             599  01-01-97  7:51p bartlett.g
LAGMAT   G             103  07-22-98  9:29p lagmat.g

Datafiles to calculate p-values (used by cluster.g):
ASYMDST0 DAT        39,600  08-03-98  5:49a asymdst0.dat
ASYMDST0 DHT           920  08-03-98  5:49a asymdst0.dht
ASYMDST1 DAT        39,600  08-03-98  5:49a asymdst1.dat
ASYMDST1 DHT           920  08-03-98  5:49a asymdst1.dht

Datafiles:
PWT      TXT        16,908  04-16-99 11:52a pwt.txt    (ASCII)
PWT      DAT        26,880  12-15-98  6:31p pwt.dat    (Gauss)
PWT      DHT         1,024  12-15-98  6:31p pwt.dht    (Gauss header file)
BERNDURL TXT         7,188  04-16-99 11:54a berndurl.txt (ASCII)
BERNDURL DHT           256  03-19-97  7:19p berndurl.dht (Gauss)
BERNDURL DAT        11,264  03-19-97  7:19p berndurl.dat (Gauss header file)

PCAG0    TXT
PCAG1    TXT
PCMAN0   TXT
PCMAN1   TXT
PCMS0    TXT
PCMS1    TXT
PCNMS0   TXT
PCNMS1   TXT


An earlier version of this file:
README   TXT                12-15-98  6:37p readme.txt

Both data files contain per capita RGDP levels in ppp dollars with base year
1980 for the Bernard and Durlauf data and 1985 for the PWT data.

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


GAUSS-procedures: Hopefully, the comments suffice when using them. They have
been tested in GAUSS 3.2.1 for DOS. The main procedure is CLUSTER.G. The
example programs illustrate how to use it.

If you have any questions or problems, please contact Bart Hobijn
(hobijnb@fasecon.econ.nyu.edu).

