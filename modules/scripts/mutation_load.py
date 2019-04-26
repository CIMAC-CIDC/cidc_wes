#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to calculate mutation load / MB
OUTPUT: txt file
"""

import os
import sys
from optparse import OptionParser

def Calculate_Nmutations(infile, outfile = "outfile.txt"):
    ofile = open(outfile, "w")
    i=0
    with open(infile,'r') as fin:
        for line in fin:
            if(line[0]!='#') or (line[0]!='Hugo_Symbol'):
                i=i+1
    if 'MDA' in infile:
        rate=float(i)/28969900
    elif 'Mocha' in infile or  'merge' in infile:
        rate=float(i)/34208715
    elif 'Broad' in infile:
        rate=float(i)/33062838
    else:
        rate=float(i)/30000000 #30000000 is general number。
    #return rate*1000000
    ofile.write("the mutation in protein coding\n" + "the mutation load/MB : " + str(rate))
    ofile.close()


def main():
    usage = "USAGE: %prog -v [maf_file.maf] -o [output file.txt]"
    optparser = OptionParser(usage=usage)                    
    optparser.add_option("-v", "--maf", help="maf file to filter")
    optparser.add_option("-o", "--out", help="output file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.maf or not options.out:
        optparser.print_help()
        sys.exit(-1)
        
    Calculate_Nmutations(options.maf, options.out)

if __name__=='__main__':
    main()
