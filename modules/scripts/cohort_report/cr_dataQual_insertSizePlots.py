#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate insert size plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser
from cr_utils import processSampleJson

def getCols(runs, attrib):
    """Returns the column values over all of the samples"""
    ls = []
    for r in runs:
        ls.append((r['tumor']['id'], r['tumor'][attrib]))
        if 'normal' in r:
            ls.append((r['normal']['id'], r['normal'][attrib]))
    return ls

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output csv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    attrs = ['gc_content', 'insert_size']
    runs = [processSampleJson(f, 'alignment', attrs) for f in options.files]

    #FOR line graphs, the hdr should be X, Sample1, Sample2,...,SampleN
    #Where the first col, X represents the x vals, and the other cols repr
    #the Y-vals from each sample (at point x)
    samples = getCols(runs, "insert_size")

    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = ['X']
    hdr.extend([s[0] for s in samples])
    out.write("%s\n" % ",".join(hdr))

    #The X axis for insert size runs from 2 to N inclusive
    #find min len
    lens = list(map(lambda s: len(s[1]), samples))
    minLen = min(lens)
    xaxis = [i+2 for i in range(minLen - 2)]
    #print(xaxis)
    for (i,xval) in enumerate(xaxis):
        row = [str(xval)]
        row.extend([str(s[1][i]) for s in samples])
        out.write("%s\n" % ",".join(row))
    out.close()

if __name__ == '__main__':
    main()
