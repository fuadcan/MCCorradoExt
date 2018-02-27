
    cls;
    print "        example for writing/reading Excel files";
    print;
    print "This program takes an existing dataset (sci.dat) on the \\EXAMPLES";
    print "directory, writes some of the data out to an Excel file, and then";
    print "reads the data back in, and displays it.";

#ifUNIX
    print;
    print "This example requires the Win/32 Operating System";
#else

    print;
    print "Press any key to continue"; call keyw;

    cls;
    gpath = sysstate(2,1);
    gname = gpath $+ "\\examples\\sci.dat";
    fname = gpath $+ "\\examples\\sci.xls";
    names = getname(gname);
    call exportf(gname,fname,names[1:6]);

    print; 
    print "Press any key to continue"; call keyw;
    {x,nm} = import(fname,0,0);
    print; 

    print "Press any key to continue"; call keyw;
    cls;
    print "Imported vectors " $nm;

    print;
    format 8,4;
    print "Data (first 10 obs)";
    print x[1:10,.];

#endif

