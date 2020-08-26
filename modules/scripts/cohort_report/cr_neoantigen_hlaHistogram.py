#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import json
import numpy as np
import pandas as pd
from optparse import OptionParser
from collections import Counter

def getHits(hla_list):
    """Given a list of HLA alleles, e.g. HLA-A alleles,
    returns a list of alleles that occur more than once"""
    count = Counter(hla_list)
    hits = [k for k, v in count.items() if v > 1]
    #cull the list
    ret = list(filter(lambda a: a in hits, hla_list))
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

    #group columns together
    cols = {"HLA-A": ('A-1', 'A-2'),
            "HLA-B": ('B-1', 'B-2'),
            "HLA-C": ('C-1', 'C-2'),
            'HLA-Dpb1': ('Dpb1-1', 'Dpb1-2'),
            'HLA-Dqb1': ('Dqb1-1', 'Dqb1-2'),
            'HLA-Drb1': ('Drb1-1', 'Drb1-2')}
    hits = [getHits(df[c1].tolist() + df[c2].tolist()) for (c1,c2) in cols.values()]
    maxRows = max(list(map(lambda ls: len(ls), hits)))
    #extend this to max rows
    for ls in hits:
        diff = maxRows - len(ls)
        if diff > 0:
            ls.extend([np.nan for i in range(diff)])
    #print(hits)

    #Create new df and write to output
    #print(list(zip(hits))[0])
    df = pd.DataFrame(list(zip(*hits)), columns=cols.keys())
    df.to_csv(options.output)

if __name__ == '__main__':
    main()
