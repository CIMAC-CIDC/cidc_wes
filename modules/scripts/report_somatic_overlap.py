#!/usr/bin/env python3
"""Len Taing 2020 (TGBTG)
Script to generate wes TMB information"""

import os
import sys
import math

from optparse import OptionParser

#MOVED to somatic_genStats.py
# def countBedBases(ffile):
#     """Count the number of base pairs in a bed file"""
#     count = 0
#     f = open(ffile)
#     for l in f:
#         tmp = l.split("\t")
#         (start, end) = (int(tmp[1]), int(tmp[2]))
#         count = count + end - start
#     f.close()
#     return count

def main():
    usage = "USAGE: %prog -f [vcf compare file from germline] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="vcf compare file")
    optparser.add_option("-r", "--run", help="run name")
    #optparser.add_option("-t", "--target", help="target region bed file")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.run or not options.output:
        optparser.print_help()
        sys.exit(-1)
    
    f = open(options.file)
    #skip first 7 lines
    for i in range(7):
        tmp = f.readline()
    tmp = f.readline().strip().split()
    normal_uniq = int(tmp[1])
    tmp = f.readline().strip().split()
    tumor_uniq = int(tmp[1])
    tmp = f.readline().strip().split()
    common = int(tmp[1])
    
    normal_tmb = common + normal_uniq
    tumor_tmb = common + tumor_uniq
    #REPORT % tumor variants that are in germline
    overlap = "%.2f" % (float(common)/tumor_tmb*100.0)
    f.close()

    out = open(options.output, 'w')
    hdr = ['Run', 'Tumor', 'Normal', 'Common', '% overlap']
    out.write("%s\n" % "\t".join(hdr))
    out.write("%s\n" % "\t".join([options.run, str(tumor_tmb), str(normal_tmb), str(common), overlap]))
    out.close()
    
if __name__ == '__main__':
    main()
