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
    optparser.add_option("-n", "--normal", help="CIMAC id of normal sample")
    optparser.add_option("-t", "--tumor", help="CIMAC id of tumor sample")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.run or not options.output or not options.normal or not options.tumor:
        optparser.print_help()
        sys.exit(-1)
    
    f = open(options.file)
    for l in f:
        #skip comments
        if l.startswith("#"):
            continue
        #VN lines are the lines we're interested in
        if l.startswith("VN"):
            #split into parts: VN, count, filepath(s), percent
            tmp = l.split()
            if len(tmp) > 4: #more than 4 parts, it's the common line
                common = int(tmp[1])
            else: #it's a sample specific line
                #parse out the CIMAC id, rely on the fact that the path is
                #in the form: analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz
                sample_name = tmp[2].split("/")[2]
                if sample_name == options.normal:
                    normal_uniq = int(tmp[1])
                else:
                    tumor_uniq = int(tmp[1])
    f.close()
    #print(common, normal_uniq, tumor_uniq)
    
    normal_variants = common + normal_uniq
    tumor_variants = common + tumor_uniq
    #REPORT % tumor variants that are in germline
    overlap = "%.2f" % (float(common)/tumor_variants*100.0)

    out = open(options.output, 'w')
    hdr = ['Run', 'Tumor', 'Normal', 'Common', '% overlap']
    out.write("%s\n" % "\t".join(hdr))
    out.write("%s\n" % "\t".join([options.run, str(tumor_variants), str(normal_variants), str(common), overlap]))
    out.close()
    
if __name__ == '__main__':
    main()
