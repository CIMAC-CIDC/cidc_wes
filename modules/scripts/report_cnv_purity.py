#!/usr/bin/env python3
"""Len Taing 2020 (TGBTG)
Script to generate report hla information"""

import os
import sys
import math

from optparse import OptionParser

def pp(s):
    """prettyprint numbers"""
    if s == 'NA':
        s = '0.0'
    return "%.4f" % float(s)

def main():
    usage = "USAGE: %prog -f [purity results file] -r [run name] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="purity result file")
    optparser.add_option("-r", "--run", help="run name")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.run or not options.output:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.file)
    hdr = f.readline().strip().split("\t")
    ls = f.readline().strip().split("\t")
    tmp = dict(zip(hdr,ls))
    f.close()
    
    out = open(options.output, 'w')
    out.write("%s\n" % "\t".join(['Run','Tumor Purity', 'Tumor Ploidy']))
    out.write("%s\n" % "\t".join([options.run, pp(tmp['purity']),
                                  pp(tmp['ploidy'])]))
    out.close()
    
if __name__ == '__main__':
    main()
