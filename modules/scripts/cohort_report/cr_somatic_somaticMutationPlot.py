#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

#Attribs to read in
#LEAVE off the plot for now
_attrs = ['mutation_summary']#, 'transition_matrix', 'tmb', 'functional_summary']

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    somatic = tmp['somatic']
    #print(somatic)
    run = {'id': tmp['id']}
    for a in _attrs:
        run[a] = somatic[a]
    #print(run)
    #sys.exit()
    return run

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def calcTiTv(trans_mat):
    """Calculates the Transition/Transversion (Ti/Tv) ratio"""
    _transitions = ["AG", "GA", "CT", "TC"]
    total = 0
    ti_ct = 0
    for (k, row) in trans_mat.items():
        for (base, n) in row.items():
            total += n
            tmp = k + base
            if tmp in _transitions:
                ti_ct += n
    #print(total, ti_ct)
    ratio = "%.3f" % (float(ti_ct)/float(total-ti_ct))
    return ratio
    
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
    runs = [processJson(f) for f in options.files]

    #build mutation summaries dict
    mutation_summaries = {}
    
    out = open(options.output, "w")
    hdr = ['Run']
    cols = list(runs[0]['mutation_summary'].keys())
    #remove total-
    cols.remove('total')
    hdr.extend(cols)
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        tmp = [r['id']]

        vals = [str(r['mutation_summary'][c]) for c in cols]
        tmp.extend(vals)
        out.write("%s\n" % ",".join(tmp))
    out.close()

if __name__ == '__main__':
    main()
