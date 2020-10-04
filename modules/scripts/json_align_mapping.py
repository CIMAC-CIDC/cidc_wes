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
    #READ next three lines for total, mapped, dedup
    total = int(f.readline().strip().split(" ")[0])
    mapped = int(f.readline().strip().split(" ")[0])
    dedup = int(f.readline().strip().split(" ")[0])
    f.close()

    ret = {'total_reads': total, 'mapped_reads': mapped, 'dedup_reads': dedup}

    #print(ret)
    f.close()
    return ret

def main():
    usage = "USAGE: %prog -r run name -t tumor file -n normal file -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-t", "--tumor", help="tumor file", default=None)
    optparser.add_option("-n", "--normal", help="tumor file", default=None)
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
    
    js_out['tumor'] = {'id': tmr_id, 'alignment': tmr}
    js_out['normal'] = {'id': nrm_id, 'alignment': nrm}

    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
