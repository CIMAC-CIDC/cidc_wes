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

#Attribs to read in
#LEAVE off the plot for now
_attrs = ['A-1', 'A-2',
          'B-1', 'B-2',
          'C-1', 'C-2',
          'DPB1-1', 'DPB1-2',
          'DQB1-1', 'DQB1-2',
          'DRB1-1', 'DRB1-2']

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
    runs = [processSampleJson(f, 'hla', _attrs) for f in options.files]

    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = map(lambda x: prettyprint(x), _attrs)
    hdr = ['Sample'] + list(hdr)
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        #Print both tumor and normal rows
        tmr = r['tumor']
        tmp = [str(tmr[a]) for a in _attrs]
        first_row = [tmr['id']]
        first_row.extend(tmp)
        out.write("%s\n" % ",".join(first_row))

        if 'normal' in r:
            nrm = r['normal']
            tmp = [str(nrm[a]) for a in _attrs]
            secnd_row = [nrm['id']]
            secnd_row.extend(tmp)
            out.write("%s\n" % ",".join(secnd_row))
    out.close()

if __name__ == '__main__':
    main()
