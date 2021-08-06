#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
import base64
from optparse import OptionParser


def main():
    usage = "USAGE: %prog -f file -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-f", "--in_file", help="input maf file", default=None)
    optparser.add_option("-j", "--tri_mtrx_file", help="input trinucleotid matrix file", default=None)
    optparser.add_option("-t", "--tmb_file", help="somatic mutation summary file which includes TMB as last col", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.in_file or not options.tri_mtrx_file or not options.tmb_file or  not options.output:
        optparser.print_help()
        sys.exit(-1)


    f = open(options.in_file)
    s = f.read()
    s_byte = s.encode('utf-8')
    #NOTE: .decode is needed to convert the bytes back to a string
    s_b64 = base64.b64encode(s_byte).decode('utf-8')
    f.close()

    #read in tri_nucleotide matrix file
    f = open(options.tri_mtrx_file)
    tri_mtrx= json.load(f)
    f.close()

    #Get the TMB from the last col of the tmb_file
    f = open(options.tmb_file)
    hdr = f.readline().strip().split(",")
    tmb = float(f.readline().strip().split(",")[-1])
    f.close()

    js_out = {'id': options.run, 'somatic': {'filtered_maf_file':"%s" % s_b64, 'tri_matrix': tri_mtrx, 'tmb': tmb}}

    out = open(options.output, 'w')
    out.write(json.dumps(js_out))
    out.close()
    
if __name__=='__main__':
    main()
