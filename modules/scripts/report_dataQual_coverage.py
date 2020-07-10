#!/usr/bin/env python3
"""Len Taing 2020 (TGBTG)
Script to generate wes coverage information"""

import os
import sys
import math

from optparse import OptionParser

def millify(n):
    """Given a large int n, returns a string representation of n in human-
    readable form
    ref: https://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
    """
    millnames = ['',' K',' M',' B',' T']

    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                    int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.1f}{}'.format(n / 10**(3 * millidx), millnames[millidx])

def main():
    usage = "USAGE: %prog -f [all sample summaries file] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="all sample summaries file")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    
    f = open(options.file)
    out = open(options.output, 'w')
    hdr = f.readline().strip().split("\t")
    #SET the hdr cols
    hdr[0] = "Sample"
    hdr[1] = "Total Reads"
    hdr[2] = "Mean Depth"
    hdr[3] = "Q1 Depth"
    hdr[4] = "Median Depth"
    hdr[5] = "Q3 Depth"
    hdr[6] = "% Bases Above 50x"
    out.write("%s\n" % "\t".join(hdr))
    for l in f:
        tmp = l.strip().split("\t")
        #Millify the total reads
        if int(tmp[1]) >= 1000:
            tmp[1] = millify(tmp[1])

        out.write("%s\n" % "\t".join(tmp))
    f.close()
    out.close()
    
if __name__ == '__main__':
    main()
