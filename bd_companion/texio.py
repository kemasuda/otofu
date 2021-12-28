#%% latex IO
def add_target(d, name, tfile, keys1, keys2, pdic):
    ncomp = min(len(d), 5) # at most 5 companions

    amp = " & "
    end = " \\\\"
    nodata = "\\nodata"

    shead = r"\subsection{%s}"%name
    thead = r"\begin{deluxetable}{lc}[h!] \caption{%s} \tablehead{\colhead{} & \colhead{target}} \startdata"%name
    ttail = r"\enddata\end{deluxetable}"

    #print (shead, file=tfile)
    print (" ", file=tfile)
    print (thead, file=tfile)
    for i,k in enumerate(keys1):
        if k=='name' or k=='spt_opt' or k=='spt_ir' or k=='exoplanet':
            val = str(d.iloc[0][k])
        else:
            val = '$%.4f$'%float(d.iloc[0][k])
        row = pdic[k] + amp + val
        #if k!='Ks_2MASS':
        if k!='z_P1':
            row += end
        print (row, file=tfile)

    print (ttail, file=tfile)

    #"""
    thead = r"\begin{deluxetable}{l@{\hspace{.5cm}}%s}[h!] \caption{Gaia sources near %s} \tablehead{\colhead{} & \multicolumn{%d}{c}{sources}} \startdata"%("c"*ncomp, name, ncomp)
    print (" ", file=tfile)
    print (thead, file=tfile)
    for i,k in enumerate(keys2):
        row = pdic[k]
        for j in range(ncomp):
            if k=='source_id':
                try:
                    val = str(int(float(d.iloc[j][k])))
                except:
                    val = "N/A"
            else:
                val = '$%.4f$'%d.iloc[j][k]
            row += amp + val
        if k!='parallax':
            row += end
        print (row, file=tfile)
    print (ttail, file=tfile)
    #"""

    fig = r"\begin{figure}[h!] \epsscale{1.2} \plotone{%s.png}\caption{%s}\end{figure}"%(name, name)
    print (" ", file=tfile)
    print (fig, file=tfile)
    print ("\\clearpage", file=tfile)
    print (" ", file=tfile)
    print (" ", file=tfile)
