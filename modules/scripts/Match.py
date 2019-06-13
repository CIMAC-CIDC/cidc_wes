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

    #Input_file = options.input
    f = open(options.input, 'r')
    lines = f.readlines()
    sum = 0
    for line in lines:
        if ("non-matching" in line) or ("main file" in line) or ("main file" in line):
            sum += int(line.split()[1])
    if sum == 0:
        print ("match")
    else:
        print("mismatch")

if __name__=='__main__':
    main()
