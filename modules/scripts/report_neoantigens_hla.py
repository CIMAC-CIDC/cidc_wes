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

#DEPRECATED
# def parseXHLA(xhla_f):
#     f = open(xhla_f)
#     xhla_out = json.load(f)
#     f.close()

#     #build classII alleleles
#     #ONLY add class II alleles--i.e. ones that start with "D"
#     classII = [a for a in xhla_out['hla']['alleles'] if a.startswith("D")]
#     #print(classII)
#     return classII

def _cleanAllele(s):
    """takes a strine HLA-DRB1*13:02:01 to DRB1*13:02"""
    #Remove the HLA-                                                        
    s = s.split("-")[1]
    #Remove the last :01                                                    
    s = ":".join(s.split(":")[:-1])
    #print(s)                                                               
    return s

def parseHLA_HD(hlahd_f, refurn_ls = False):
    """PARSES the results/{name}_final.txt output from hla-hd"""
    #PARSE hlahd txt file...
    c2_alleles = {}
    #NOTE read only A,B,C,DRB1,DQA1,DQB1,DPA1,DPB1 (first 8 lines)
    #REST: DMA,DMB,DOA,DOB,DRA,DRB2,DRB3,DRB4,DRB5,DRB6,DRB7,DRB8,DRB9,
    #DPA2,E,???,G,H,J,K,L,???,V,???,Y
    if os.path.exists(hlahd_f):
        f = open(hlahd_f)
        for i in range(8): #first 8 lies
            l = f.readline().strip()
            if i > 2 and not l.startswith('Could'):#ksip 'Couldn't read result file.' and A,B,C alleles
                tmp = l.split("\t")
                #Coerce HLA-DRB1*13:02:01 to DRB1*13:02
                c2_alleles[tmp[0]] = [_cleanAllele(a) for a in tmp[1:3] if a != "-"]
        f.close()
    #print(c2_alleles)

    if refurn_ls:
        #RETURN the alleles as a list
        c2_list = [",".join(a) for a in c2_alleles.values()]
        #print(c2_list)
        return c2_list
    else:
        return c2_alleles

def main():
    usage = "USAGE: %prog -n [optitype/hlahd result file for normal] -t [optitype/hlahd result file for tumor] -s [sample names, normal first] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-n", "--normal_opti", help="hla result file normal")
    optparser.add_option("-m", "--normal_hlahd", help="hla result file normal")
    optparser.add_option("-t", "--tumor_opti", help="hla result file tumor")
    optparser.add_option("-u", "--tumor_hlahd", help="hla result file tumor")
    optparser.add_option("-s", "--names", help="comma separated sample names, normal first")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.tumor_opti or not options.names or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #normal_files = options.normal.split(",") if ',' in options.normal else [options.normal]
    #tumor_files = options.tumor.split(",") if ',' in options.tumor else [options.tumor]
    samples = options.names.split(",")

    normal_classI = parseOptitype(options.normal_opti) if options.normal_opti else None
    normal_classII = parseHLA_HD(options.normal_hlahd, True) if options.normal_hlahd else None

    tumor_classI = parseOptitype(options.tumor_opti)
    tumor_classII = parseHLA_HD(options.tumor_hlahd, True) if options.tumor_hlahd else None

    out = open(options.output,"w")
    hdr = ["Sample", "A1", "A2", "B1", "B2", "C1", "C2"]
    out.write("%s\n" % "\t".join(hdr))

    #Check if this is tumor-only--assume that if normal is None, then tmr only
    if not options.normal_opti:
        tumor_classI.insert(0, samples[0])
        out.write("%s\n" % "\t".join(tumor_classI))
        if tumor_classII:
            tumor_classII.insert(0, '&nbsp;')
            #ALSO pad the end
            tumor_classII.append('&nbsp;')
            out.write("%s\n" % "\t".join(tumor_classII))
        out.close()
    else:
        normal_classI.insert(0, samples[0])
        out.write("%s\n" % "\t".join(normal_classI))
        if normal_classII:
            normal_classII.insert(0, '&nbsp;')
            normal_classII.append('&nbsp;')
            out.write("%s\n" % "\t".join(normal_classII))

        tumor_classI.insert(0, samples[1])
        out.write("%s\n" % "\t".join(tumor_classI))
        if tumor_classII:
            tumor_classII.insert(0, '&nbsp;')
            tumor_classII.append('&nbsp;')
            out.write("%s\n" % "\t".join(tumor_classII))
        out.close()
    
if __name__ == '__main__':
    main()
