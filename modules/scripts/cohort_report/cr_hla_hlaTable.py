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
_attrs = ['A-1', 'A-2',
          'B-1', 'B-2',
          'C-1', 'C-2',
          'DPB1-1', 'DPB1-2',
          'DQB1-1', 'DQB1-2',
          'DRB1-1', 'DRB1-2']

def getSampleInfo(d, isTumor):
    """Given a dictionary of json records, returns a record of tumor OR
    normal sample"""
    sample = d['tumor'] if isTumor else d['normal']
    hla = sample['hla']
    
    tmp = {'id': sample['id']}
    for a in _attrs:
        tmp[a] = hla.get(a)
    return tmp

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    tumor = getSampleInfo(tmp, True)
    normal = getSampleInfo(tmp, False)
    run = {'id': tmp['id'], 'tumor':tumor, 'normal': normal}
    #print(run)
    return run

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

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

        nrm = r['normal']
        tmp = [str(nrm[a]) for a in _attrs]
        secnd_row = [nrm['id']]
        secnd_row.extend(tmp)
        #print(first_row)
        #print(secnd_row)
        out.write("%s\n" % ",".join(first_row))
        out.write("%s\n" % ",".join(secnd_row))
    out.close()

if __name__ == '__main__':
    main()
