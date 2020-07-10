#!/usr/bin/env python3

"""
script to calculate clonality--from jingxin fu
INPUT: pyclone cluster table, e.g. {sample}_table.tsv
OUTPUT: clonality value
"""

import os
import sys
from optparse import OptionParser


import pandas as pd
import numpy as np

def estClonality(df):
    P = df['size'] / df['size'].sum()
    n = df.shape[0]
    return 1 + (P*np.log2(P)).sum()/np.log2(n)

def calcClonality(table_file):
    """Simple wrapper to call estClonality, but handles files not data frames
    This is so I can make a call from within report_level2
    """
    pd_table = pd.read_csv(table_file, delimiter="\t", encoding="utf-8")
    clonality_val = round(estClonality(pd_table), 3)
    return clonality_val

def main():
    usage = "USAGE: %prog -f [pyclone cluster table.tsv file]"
    optparser = OptionParser(usage=usage)                    
    optparser.add_option("-f", "--file", help="pyclone cluster table.tsv file")
    optparser.add_option("-r", "--run", help="run name")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.run or not options.output:
        optparser.print_help()
        sys.exit(-1)

    clonality = str(calcClonality(options.file))

    out = open(options.output, 'w')
    out.write("%s\n" % "\t".join(['Run','Tumor Clonality']))
    out.write("%s\n" % "\t".join([options.run,clonality]))
    out.close()

if __name__=='__main__':
    main()
