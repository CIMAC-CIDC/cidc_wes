#!/usr/bin/env python 
"""
This Script is for judgment of germline matching situation of samples, 
The input parameter is the log file which is generated from vcftools diff.
"""
import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -i [_diff.log]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-i", "--input", help="path to _diff.log file")
    (options, args) = optparser.parse_args(sys.argv)
    if not options.input :
        optparser.print_help()
        sys.exit(-1)

    f = open(options.input, 'r')
    for line in f:
        if ("non-matching" in line): 
            if int(line.split(" ")[1]) == 0:
                print("match")
            else:
                print("mismatch")
            sys.exit() #break out--we're done

if __name__=='__main__':
    main()
