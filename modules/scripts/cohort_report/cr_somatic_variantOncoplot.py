#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate binary matrix of top shared somatic variants"""

import os
import sys
import json
import numpy as np
import pandas as pd
import base64
from optparse import OptionParser
from collections import Counter

_attrs = ['filtered_maf_file']

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    somatic = tmp['somatic']
    #print(somatic)
    run = {'id': tmp['id']}
    for a in _attrs:
        run[a] = somatic[a]
    #print(run)
    #sys.exit()
    return run

def processMaf(maf_f, run, hits_d):
    """given a maf file, a run name, and a dictionary to modify-
    modifies the hits_d so that each key = HGVSp and 
    value = list of runs that contain that variant
    """
    #skip first two lines- version and hdr
    hdr = maf_f[1].split("\t")
    for l in maf_f[2:]:
        tmp = dict(zip(hdr, l.split("\t")))
        hgvsp = tmp.get('HGVSp', None)
        hugo = tmp.get('Hugo_Symbol', None)
        if hgvsp and hugo:
            variant ="%s:%s" % (hugo, hgvsp)
            if variant in hits_d: 
                #check to ensure run is added only once
                if run not in hits_d[variant]:
                    hits_d[variant].append(run)
            else:
                #new variant
                hits_d[variant] = [run]

def getHitsList(run, variantsLs):
    """Given a run name and a list of variants [(HGVSp, [...list of runs])]
    returns a binary list --1 if the run has the variant, 0 otherwise
    """
    tmp = list(map(lambda x: "1" if run in x[1] else "0", variantsLs))
    return tmp

def main():
    usage = "USAGE: %prog -f [wes json file1] -f [wes json file2]... -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="wes json file")
    optparser.add_option("-o", "--output", help="output csv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #variantsHits is a dictionary of variants based on the HGVSp;
    #key = HGVSp, value = list of runs that contain that variant
    variantsHits = {}
    #READ in maf files and build up binary matrix
    for (run, maf_b64) in [(r['id'], r['filtered_maf_file']) for r in runs]:
        maf = base64.b64decode(maf_b64).decode("utf-8").split("\n")
        #process maf file and build up variantsHits
        processMaf(maf, run, variantsHits)

    sorted_variants = sorted(variantsHits.items(), key=lambda x: len(x[1]), reverse=True)

    #TAKE top 100 variants
    sorted_variants = sorted_variants[:100]

    out = open(options.output, "w")
    hdr = ["Run"]
    hdr.extend(list(map(lambda x: x[0], sorted_variants)))
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        run_name = r['id']
        hits = getHitsList(run_name, sorted_variants)
        row = [run_name]
        row.extend(hits)
        out.write("%s\n" % ",".join(row))
    out.close()

if __name__ == '__main__':
    main()
