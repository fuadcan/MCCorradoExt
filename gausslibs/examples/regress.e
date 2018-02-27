new;
cls;

print "The regress.e example illustrates how a keyword procedure";
print "can be implemented.  After you have run regress.e which";
print "compiled code to enable the keyword named regress.";
print "Follow the on screen instructions.";
print;
print "Press any key to continue"; call keyw;
cls;

where = cdir(0);
dir = sysstate(2,0) $+ "examples";
chdir ^dir;

keyword regress(s);

    local dep,ind,tok,dataset;
    dep = 0;
    ind = 0;

    { dep,s } = token(s);
    { tok,s } = token(s);

    if (upper(tok) $/= "ON");
        goto usage;
    endif;

    do until (upper(tok) $== "FROM");
        { tok,s } = token(s);
        if s $== "";
            goto usage;
        endif;
        ind = ind|tok;
    endo;

    ind = trimr(ind,1,1);
    dataset = s;
    call ols(dataset,dep,ind);
    retp;
Usage:

    print;
    print "***************************************************************";
    print "*                                                             *";
    print "*  Usage:                                                     *";
    print "*                                                             *";
    print "*    REGRESS <depvar> ON <indvar list> FROM <dataset>;        *";
    print "*                                                             *";
    print "*  Example:                                                   *";
    print "*                                                             *";
    print "*    regress pub6 on pub3 job mmale from sci;                 *";
    print "*                                                             *";
    print "***************************************************************";

endp;

__miss = 1;
regress pub6 on pub3 job mmale from sci;

chdir ^where;

print;
print;

print "Now type regress at the prompt, and hit Enter. The usage will be";
print "explained.  You also need to be in the examples directory, so ";
print "that the GAUSS data set used in this example is in the same";
print "directory";
print;
print "Your current location is ";;
print cdir(0);

