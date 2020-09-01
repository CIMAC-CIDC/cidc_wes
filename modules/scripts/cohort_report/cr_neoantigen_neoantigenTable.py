#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

#Attribs to read in
_attrs = ['Gene_Name', 'Mutation', 'Protein_Position', 'HGVSc', 'HGVSp', 'HLA_Allele', 'MT_Epitope_Seq', 'MT_IC50', 'WT_IC50', 'Fold_Change', 'Tumor_DNA_Depth', 'Tumor_DNA_VAF', 'Score']


def roundFloat(entry):
    """Tries to round the float fields to two decimal places"""
    _float_attrs = ['MT_IC50','WT_IC50','Fold_Change','Tumor_DNA_VAF','Score']

    #print(entry)
    for a in _float_attrs:
        if a in entry:
            entry[a] = "%.2f" % float(entry[a])
    return entry

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()


    ls = [roundFloat(e) for e in tmp['neoantigen']]
    run = {'id': tmp['id'], 'neoantigen': ls}
    #print(run)
    return run

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output tsv file")
    optparser.add_option("-j", "--json", help="json output file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #NOTE: for output, just place the top hit; for json dump the entire list

    neoantigens = {}
    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = _attrs
    hdr = ['Run'] + list(hdr)
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        #Print both tumor and normal rows
        top_hit = [str(r['neoantigen'][0][a]) for a in _attrs]
        row = [r['id']]
        row.extend(top_hit)
        #print(row)

        out.write("%s\n" % ",".join(row))
        #build the json-
        neoantigens[r['id']] = r['neoantigen']
    out.close()

    json_out = open(options.json, "w")
    json_out.write(json.dumps(neoantigens))
    json_out.close()


if __name__ == '__main__':
    main()
