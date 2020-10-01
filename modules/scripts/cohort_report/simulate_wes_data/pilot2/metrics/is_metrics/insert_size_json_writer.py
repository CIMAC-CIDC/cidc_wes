#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser


def parseFile(gc_file):
    """Reads in a _gc_metrics.txt file- 4th col"""
    f = open(gc_file)
    #Burn first 5 lines
    for i in range(5):
        hdr = f.readline()

    ls = []
    for l in f:
        if l.strip():
            tmp = l.strip().split("\t")
            ls.append(int(tmp[1])) #grab 2nd col
    f.close()
    return ls

def main():
    usage = "USAGE: %prog -r run_name -t tumor_file -n normal_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-t", "--tumor", help="tumor file", default=None)
    optparser.add_option("-n", "--normal", help="normal file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.tumor or not options.normal or not options.output:
        optparser.print_help()
        sys.exit(-1)

    js_out = {'id': options.run}
    tmr = parseFile(options.tumor)
    #GET tumor sample name
    fname = options.tumor.split("/")[-1]
    tmr_id = fname.split("_")[0]
    
    nrm = parseFile(options.normal)
    #GET normal sample name
    fname = options.normal.split("/")[-1]
    nrm_id = fname.split("_")[0]
    
    js_out['tumor'] = {'id': tmr_id, 'alignment': {'insert_size': tmr}}
    js_out['normal'] = {'id': nrm_id, 'alignment': {'insert_size': nrm}}

    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
