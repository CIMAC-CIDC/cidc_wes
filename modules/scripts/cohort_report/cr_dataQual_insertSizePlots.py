#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate insert size plots information"""

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
    attrs = ['gc_content', 'insert_size']
    
    tmp = {'id': sample['id']}
    for a in attrs:
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

def millify(n):
    """Given a large int n, returns a string representation of n in human-
    readable form
    ref: https://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
    """
    millnames = ['',' K',' M',' B',' T']

    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                    int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.1f}{}'.format(n / 10**(3 * millidx), millnames[millidx])

def getCols(runs, attrib):
    """Returns the column values over all of the samples"""
    ls = []
    for r in runs:
        ls.append((r['tumor']['id'], r['tumor'][attrib]))
        ls.append((r['normal']['id'], r['normal'][attrib]))
    return ls

def printCell(sample, stats, attrib):
    n = sample[attrib]
    mean, std = stats[attrib]
    if n < (mean - 1.5*std): #1.5 ~ bottom 10%, 1 std ~ bottom 33%
        s = "red:%s" % millify(n)
    else:
        s = millify(n)
    return s

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
    runs = [processJson(f) for f in options.files]

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
