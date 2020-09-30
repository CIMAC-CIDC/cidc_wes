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
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-t", "--tumor_opti", help="tumor optitype results file", default=None)
    optparser.add_option("-s", "--tumor_xhla", help="tumor xhla results file", default=None)
    optparser.add_option("-n", "--normal_opti", help="normal optitype results file", default=None)
    optparser.add_option("-m", "--normal_xhla", help="normal xhla results file", default=None)

    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.tumor_opti or not options.tumor_xhla or not options.normal_opti or not options.normal_xhla or not options.output:
        optparser.print_help()
        sys.exit(-1)

    js_out = {'id': options.run}
    tmr = parseFile(options.tumor_opti, options.tumor_xhla)
    #GET tumor sample name
    fname = options.tumor_opti.split("/")[-1]
    tmr_id = fname.split("_")[0]

    nrm = parseFile(options.normal_opti, options.normal_xhla)
    #GET normal sample name
    fname = options.normal_opti.split("/")[-1]
    nrm_id = fname.split("_")[0]

    #print(total, mapped, dedup)
    js_out['tumor'] = {'id': tmr_id, 'hla': tmr}
    js_out['normal'] = {'id': nrm_id, 'hla': nrm}
    
    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
