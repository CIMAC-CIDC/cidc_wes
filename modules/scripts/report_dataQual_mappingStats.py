#!/usr/bin/env python3
"""Script to generate wes run information"""

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
    usage = "USAGE: %prog -f [mapping stats .csv file] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    
    f = open(options.file)
    out = open(options.output, "w")
    hdr = f.readline().strip().split(",")
    out.write("%s\n" % "\t".join(hdr))
    for l in f:
        tmp = l.strip().split(",")
        #Millify numbers
        tmp[1] = millify(tmp[1])
        tmp[2] = millify(tmp[2])
        out.write("%s\n" % "\t".join(tmp))
    f.close()
    out.close()
if __name__ == '__main__':
    main()
