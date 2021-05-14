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


def getCol_Stats(runs, attrib):
    """Returns the column values over all of the samples"""
    ls = []
    for r in runs:
        ls.append((r['tumor']['id'], r['tumor'][attrib]))
        if 'normal' in r:
            ls.append((r['normal']['id'], r['normal'][attrib]))
    vals = [x[1] for x in ls]
    mean = sum(vals)/len(vals)
    std = np.std(vals)
    return (ls, mean, std)

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
    attrs = ['total_reads', 'mapped_reads', 'dedup_reads', 'mean_quality_score']
    runs = [processSampleJson(f, 'alignment', attrs) for f in options.files]

    #get the col values and calculate mean and std-dev for each attr
    stats = {}
    #attrs = ['total_reads', 'mapped_reads', 'dedup_reads', 'mean_quality_score']
    attrs = ['total_reads', 'mapped_reads', 'dedup_reads']
    for a in attrs:
        stats[a] = getCol_Stats(runs, a)

    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = map(lambda x: prettyprint(x), attrs)
    hdr = ['Sample'] + list(hdr)
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        #Print both tumor and normal rows
        tmr = r['tumor']
        #tmp = [millify(tmr.get(a)) for a in attrs]
        tmp = [str(tmr[a]) for a in attrs]
        #BASE: dedup reads--the mapped and total reads are deltas
        #tmp = [str(tmr['total_reads'] - tmr['mapped_reads']),
        #       str(tmr['mapped_reads'] - tmr['dedup_reads']),
        #       str(tmr['dedup_reads'])]
        first_row = [tmr['id']]
        first_row.extend(tmp)
        out.write("%s\n" % ",".join(first_row))

        if 'normal' in r:
            nrm = r['normal']
            #tmp = [millify(nrm.get(a)) for a in attrs]
            tmp = [str(nrm[a]) for a in attrs]
            #tmp = [str(nrm['total_reads'] - nrm['mapped_reads']),
            #       str(nrm['mapped_reads'] - nrm['dedup_reads']),
            #       str(nrm['dedup_reads'])]
            secnd_row = [nrm['id']]
            secnd_row.extend(tmp)
            out.write("%s\n" % ",".join(secnd_row))
    out.close()

if __name__ == '__main__':
    main()
