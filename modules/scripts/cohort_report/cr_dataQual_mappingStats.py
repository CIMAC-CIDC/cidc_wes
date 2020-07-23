#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

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

def getStats(runs, attrib):
    ls = []
    for r in runs:
        ls.append(r['tumor'][attrib])
        ls.append(r['normal'][attrib])
    mean = sum(ls)/len(ls)
    std = np.std(ls)
    return (mean, std)

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
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #calculate mean and std-dev for each attr
    stats = {}
    attrs = ['total_reads', 'mapped_reads', 'dedup_reads', 'mean_quality_score']
    for a in attrs:
        stats[a] = getStats(runs, a)
        
    #print out data
    out = open(options.output, "w")
    hdr = ['Run','Sample','Total Reads','Mapped Reads','Dedup Reads','Mean Quality Score']
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        #Print both tumor and normal rows
        tmr = r['tumor']
        #tmp = [millify(tmr.get(a)) for a in attrs]
        tmp = [printCell(tmr, stats, a) for a in attrs]
        first_row = [r['id'],tmr['id']]
        first_row.extend(tmp)
        
        nrm = r['normal']
        #tmp = [millify(nrm.get(a)) for a in attrs]
        tmp = [printCell(nrm, stats, a) for a in attrs]
        secnd_row = ['&nbsp;',nrm['id']]
        secnd_row.extend(tmp)
        #print(first_row)
        #print(secnd_row)
        out.write("%s\n" % ",".join(first_row))
        out.write("%s\n" % ",".join(secnd_row))
    out.close()
if __name__ == '__main__':
    main()
