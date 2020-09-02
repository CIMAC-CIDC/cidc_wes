#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def getSampleInfo(d, isTumor):
    """Given a dictionary of json records, returns a record of tumor OR
    normal sample"""
    sample = d['tumor'] if isTumor else d['normal']
    alignment = sample['alignment']
    attrs = ['total_reads', 'mapped_reads', 'dedup_reads', 'mean_quality_score']
    tmp = {'id': sample['id']}
    for a in attrs:
        #millify at printout
        #tmp[a] = millify(alignment.get(a))
        tmp[a] = alignment.get(a)
    return tmp

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    tumor = getSampleInfo(tmp, True)
    normal = getSampleInfo(tmp, False)
    run = {'id': tmp['id'], 'tumor': tumor, 'normal': normal}
    return run

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

    attrs = ['mean_quality_score']

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

        nrm = r['normal']
        tmp = [str(nrm[a]) for a in attrs]
        secnd_row = [nrm['id']]
        secnd_row.extend(tmp)
        #print(first_row)
        #print(secnd_row)
        out.write("%s\n" % ",".join(first_row))
        out.write("%s\n" % ",".join(secnd_row))
    out.close()

if __name__ == '__main__':
    main()
