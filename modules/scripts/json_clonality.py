#!/usr/bin/env python
"""
2022-04-23 UPDATE: Dropping the clonality measure from pyclone {sample}_table.tsv b/c pyclone6's output changed significantly.  Instead incorporating CCF measures, pyclone6_input_file --for multisample analysis OFFLINE, and pyclone6_results_file ...just in case this is useful-

Example:
{"id": "<run name>", "copy_number": {"clonality": {'cellular_prevalences': [0.4017, 0.2348, 0.1201], 'pyclone6_input_file': '#b64 enc str', 'pyclone6_results_file': '#b64 enc str'}}}


OBSOLETE:
script to calculate clonality--from jingxin fu
INPUT: pyclone cluster table, e.g. {sample}_table.tsv
OUTPUT: clonality value
"""

import os
import sys
import subprocess
import json
import base64

from optparse import OptionParser

import pandas as pd
import numpy as np

#OBSOLETE
# def estClonality(df):
#     P = df['size'] / df['size'].sum()
#     n = df.shape[0]
#     return 1 + (P*np.log2(P)).sum()/np.log2(n)

# def calcClonality(table_file):
#     """Simple wrapper to call estClonality, but handles files not data frames
#     This is so I can make a call from within report_level2
#     """
#     pd_table = pd.read_csv(table_file, delimiter="\t", encoding="utf-8")
#     clonality_val = round(estClonality(pd_table), 3)
#     return clonality_val

def getCellPrevalence(pyclone6_results_f):
    """Given a pyclone6 results file, will return a list of the cellular 
    prevalence values for each cluster
    INPUT: **pyclone6 summary results file i.e. 
    analysis/clonality/{run}/{run}_pyclone6.results.summary.tsv
    
    Each row in the summary file represents a cluster; and they are ordered
    0, 1, ..., N
    
    RETURN: a list of cellular prevalances, where index 0 corresponds to 
    cluster 0's prevalence, index 1 to cluster1, etc.
    """
    with open(pyclone6_results_f) as f:
        hdr = f.readline()
        ls = [l.strip().split("\t")[1] for l in f]
    return ls

def base64encode(in_file):
    """encodes the contends of a file and returns a b64 string"""
    #ENCODE the input and results files
    f = open(in_file)
    s = f.read()
    s_byte = s.encode('utf-8')
    #NOTE: .decode is needed to convert the bytes back to a string
    s_b64 = base64.b64encode(s_byte).decode('utf-8')
    f.close()
    
    return s_b64
    

def main():
    usage = "USAGE: %prog -r run_name -f file -o output_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-i", "--in_file", help="pyclone input file to stash", default=None)
    optparser.add_option("-j", "--results_file", help="pyclone results file to stash", default=None)
    optparser.add_option("-k", "--summary_results_file", help="summary of pyclone results", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.in_file or not options.summary_results_file or not options.results_file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    cell_prevalences = getCellPrevalence(options.summary_results_file)
    pyclone6_input = base64encode(options.in_file)
    pyclone6_results = base64encode(options.results_file)
    
    js_out = {'id': options.run, 'copy_number': {'clonality': {'pyclone6_input_file': pyclone6_input, 'pyclone6_results_file': pyclone6_results, 'cellular_prevalences': cell_prevalences}}}
    
    out = open(options.output, 'w')
    out.write(json.dumps(js_out))
    out.close()
    
if __name__=='__main__':
    main()
