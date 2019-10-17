#!/usr/bin/env python
"""
removes anything that is a non-cannonical chrom, i.e. 1-22, X, Y, M
"""

import os
import sys
from optparse import OptionParser

            
def main():
    usage = "USAGE: %prog -f [tab delim txt] -o [output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--infile", help="file to filter")
    optparser.add_option("-o", "--out", help="output file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.infile or not options.out:
        optparser.print_help()
        sys.exit(-1)
    
    f = open(options.infile)
    out = open(options.out, "w")

    #rest of the file
    for l in f:
        tmp = l.strip().split("\t")
        if tmp[1].startswith("chr"):
            tmp[1] = "hs%s" % tmp[1][3:]
        out.write("%s\n" % "\t".join([tmp[1],tmp[2],tmp[3],tmp[5]]))
    f.close()
    out.close()

if __name__=='__main__':
    main()
