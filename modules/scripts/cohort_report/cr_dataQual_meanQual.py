#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser
from cr_utils import processSampleJson, prettyprint

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    attrs = ['mean_quality_score']
    runs = [processSampleJson(f, 'alignment', attrs) for f in options.files]


    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = map(lambda x: prettyprint(x), attrs)
    hdr = ['Sample'] + list(hdr)
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        #Print both tumor and normal rows
        tmr = r['tumor']
        tmp = [str(tmr[a]) for a in attrs]
        first_row = [tmr['id']]
        first_row.extend(tmp)
        out.write("%s\n" % ",".join(first_row))

        if 'normal' in r:
            nrm = r['normal']
            tmp = [str(nrm[a]) for a in attrs]
            secnd_row = [nrm['id']]
            secnd_row.extend(tmp)
            out.write("%s\n" % ",".join(secnd_row))
    out.close()

if __name__ == '__main__':
    main()
