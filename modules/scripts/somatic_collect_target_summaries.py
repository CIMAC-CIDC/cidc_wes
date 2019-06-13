#!/usr/bin/env python
"""Script to collect the sample summaries of each sample's _target_metrics.txt.sample_summary file

OUTPUTS to stdout:
HEADER: 
sample_id	total	mean	granular_Q1	granular_median	granular_Q3   %_bases_above_50

AND then the second row for each of the files given
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -f [FPKM FILE_2] ...-f [FPKM FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of _target_metrics.txt.sample_summary files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    #READ in the first file:
    first = options.files[0]
    f = open(first)
    hdr = f.readline().strip()
    #PRINT out header only for first file
    print(hdr)
    #print out second line
    print(f.readline().strip())
    f.close()

    #now just report the 2nd line for each subsequent file
    for ffile in options.files[1:]:
        f = open(ffile)
        hdr = f.readline().strip()
        tmp = f.readline().strip()
        print(tmp)
        f.close()

if __name__=='__main__':
    main()


