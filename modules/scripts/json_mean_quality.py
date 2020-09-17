#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser


def parseFile(in_file):
    """Reads in a mq_metrics.txt file"""
    f = open(in_file)
    #Burn first two lines
    hdr = f.readline().strip().split()
    hdr = f.readline().strip().split()
    mq = []
    for l in f:
        if l.strip():
            mq.append(float(l.strip().split("\t")[1]))
    f.close()

    ret = {'mean_quality_score': sum(mq)/len(mq)} #return mean}
    #print(ret)
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
    js_out = {'id': sample_id, 'alignment': contents}

    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
