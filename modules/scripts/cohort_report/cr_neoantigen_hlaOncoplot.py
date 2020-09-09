#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate binary matrix of top HLA alleles (for oncoplot fn)"""

import os
import sys
import json
import numpy as np
import pandas as pd
from optparse import OptionParser
from collections import Counter

def getHits(hla_list):
    """Given a list of HLA alleles, e.g. HLA-A alleles,
    returns a list of alleles and the number of times they occur
    """
    #print(hla_list)
    count = Counter(hla_list)
    #print(count)
    hits = [(k,v) for k, v in count.items()]
    #cull the list
    #print(hits)
    ret = sorted(hits, key=lambda x: x[1], reverse=True)
    #print(ret)
    #ret = list(filter(lambda a: a in hits, hla_list))
    return ret

def main():
    usage = "USAGE: %prog -f [01_HLA_table.mqc] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="01_HLA_table.mqc file")
    optparser.add_option("-o", "--output", help="output csv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    df = pd.read_csv(options.file, index_col=0)

    #join all HLAs together
    all_alleles = []
    for c in df.columns:
        all_alleles.extend(df[c].tolist())
    #print(all_alleles)
    #hits = [getHits(df[c1].tolist() + df[c2].tolist()) for (c1,c2) in cols.values()]
    #print(df.columns)
    #topHits = getHits(all_alleles)
    #print(topHits)
    #Return top 25 hits!
    topHits = list(map(lambda x: x[0], getHits(all_alleles)))[:25]
    out = open(options.output, "w")
    hdr = ["Sample"]
    hdr.extend(topHits)
    out.write("%s\n" % ",".join(hdr))
    for index, row in df.iterrows():
        alleles = [row[c] for c in df.columns]
        tmp = ["1.0" if c in alleles else "0.0" for c in topHits ]
        r = [index]
        r.extend(tmp)
        out.write("%s\n" % ",".join(r))
    out.close()

if __name__ == '__main__':
    main()
