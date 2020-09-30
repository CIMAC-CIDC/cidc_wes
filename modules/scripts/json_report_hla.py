#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser

#from report_neoantigens_hla.py
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

#from report_neoantigens_hla.py
def parseXHLA(xhla_f):
    f = open(xhla_f)
    xhla_out = json.load(f)
    f.close()

    #build classII alleleles
    #ONLY add class II alleles--i.e. ones that start with "D"
    classII = [a for a in xhla_out['hla']['alleles'] if a.startswith("D")]
    #print(classII)
    return classII


def parseFile(optitype_file, xhla_file):
    """Reads in a _gc_metrics.txt file- 4th col"""
    classI_alleles = ["A-1", "A-2", "B-1", "B-2", "C-1", "C-2"]
    classII_alleles = ["DPB1-1", "DPB1-2", "DQB1-1","DQB1-2","DRB1-1","DRB1-2"]
    classI = zip(classI_alleles, parseOptitype(optitype_file))
    classII = zip(classII_alleles, parseXHLA(xhla_file))

    #compose the two classes together
    hla = list(classI)
    hla.extend(list(classII))

    ret = dict(hla)
    #print(ret)
    return ret

def main():
    usage = "USAGE: %prog -r run_name -t tumor_file -n normal_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-p", "--optitype", help="optitype results file", default=None)
    optparser.add_option("-x", "--xhla", help="xhla results file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.optitype or not options.xhla  or not options.output:
        optparser.print_help()
        sys.exit(-1)

    hla = parseFile(options.optitype, options.xhla)
    #GET sample name
    fname = options.optitype.split("/")[-1]
    sample_id = fname.split("_")[0]

    #print(total, mapped, dedup)
    js_out = {'id': sample_id, 'hla': hla}
    
    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
