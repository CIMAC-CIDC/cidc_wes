#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser


def parseFile(in_file):
    """Reads in a coverage.txt file- 2nd line"""
    f = open(in_file)
    hdr = f.readline().strip().split()
    l = f.readline().strip().split()
    tmp = dict(zip(hdr, l))

    ret = {'total_reads': tmp['total'], 'mean_depth': tmp['mean'],
           'median_depth': tmp['granular_median'],
           'q1_depth': tmp['granular_Q1'],
           'q3_depth': tmp['granular_Q3'],
           'percent_bases_gt_50': tmp['%_bases_above_50']}

    #print(ret)
    f.close()
    return ret

def main():
    usage = "USAGE: %prog -f file -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--in_file", help="gc file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.in_file or not options.output:
        optparser.print_help()
        sys.exit(-1)


    contents = parseFile(options.in_file)
    #GET tumor sample name
    fname = options.in_file.split("/")[-1]
    sample_id = fname.split("_")[0]
    js_out = {'id': sample_id, 'coverage': contents}

    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
