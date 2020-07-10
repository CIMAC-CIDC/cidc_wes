#!/usr/bin/env python3
"""Len Taing 2020 (TGBTG)
Script to generate report hla information"""

import os
import sys
import math

import json
from optparse import OptionParser

def parseOptitype(optitype_f):
    """NOTE: cureently the optitype results.tsv looks somthing like this:
    	A1	A2	B1	B2	C1	C2	Reads	Objective
    0					C*06:04	C*06:04	4.0	3.99
    **So were' going to parse cols 1-6 and return that"""
    f = open(optitype_f)
    hdr = f.readline().strip().split("\t") #ignore for now
    classI = f.readline().strip().split("\t")[1:7] #want first 6 cols
    #print(classI)
    f.close()
    return classI

def parseXHLA(xhla_f):
    f = open(xhla_f)
    xhla_out = json.load(f)
    f.close()

    #build classII alleleles
    #ONLY add class II alleles--i.e. ones that start with "D"
    classII = [a for a in xhla_out['hla']['alleles'] if a.startswith("D")]
    #print(classII)
    return classII

def main():
    usage = "USAGE: %prog -n [optitype/xhla result file for normal] -t [optitype/xhla result file for tumor] -s [sample names, normal first] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-n", "--normal", help="hla result file normal")
    optparser.add_option("-t", "--tumor", help="hla result file tumor")
    optparser.add_option("-s", "--names", help="comma separated sample names, normal first")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.normal or not options.tumor or not options.names or not options.output:
        optparser.print_help()
        sys.exit(-1)

    normal_files = options.normal.split(",") if ',' in options.normal else [options.normal]
    tumor_files = options.tumor.split(",") if ',' in options.tumor else [options.tumor]
    samples = options.names.split(",")

    normal_classI = parseOptitype(normal_files[0])
    normal_classII = parseXHLA(normal_files[1]) if len(normal_files) > 1 else None

    tumor_classI = parseOptitype(tumor_files[0])
    tumor_classII = parseXHLA(tumor_files[1]) if len(tumor_files) > 1 else None

    out = open(options.output,"w")
    hdr = ["Sample", "A1", "A2", "B1", "B2", "C1", "C2"]
    out.write("%s\n" % "\t".join(hdr))
    normal_classI.insert(0, samples[0])
    out.write("%s\n" % "\t".join(normal_classI))
    if normal_classII:
        normal_classII.insert(0, '&nbsp;')
        out.write("%s\n" % "\t".join(normal_classII))

    tumor_classI.insert(0, samples[1])
    out.write("%s\n" % "\t".join(tumor_classI))
    if tumor_classII:
        tumor_classII.insert(0, '&nbsp;')
        out.write("%s\n" % "\t".join(tumor_classII))
    out.close()
    
if __name__ == '__main__':
    main()
