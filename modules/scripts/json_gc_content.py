#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser


def parseGCContent(gc_file):
    """Reads in a _gc_metrics.txt file- 4th col"""
    f = open(gc_file)
    #Burn first two lines
    hdr = f.readline()
    hdr = f.readline()
    ls = []
    for l in f:
        if l.strip():
            tmp = l.strip().split("\t")
            ls.append(tmp[3]) #grab 4th col
    f.close()
    return ls

def main():
    usage = "USAGE: %prog -r run_name -t tumor_file -n normal_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--gc_file", help="gc file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.gc_file or not options.output:
        optparser.print_help()
        sys.exit(-1)


    gc = parseGCContent(options.gc_file)
    #GET tumor sample name
    fname = options.gc_file.split("/")[-1]
    sample_id = fname.split("_")[0]
    js_out = {'id': sample_id, 'alignment': {'gc_content': gc}}

    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
