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
#LEAVE off the plot for now
_attrs = ['clonality', 'purity', 'ploidy', 'dipLogR']#, 'cnv_plot_file', 'cnv_file']

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    cnv = tmp['copy_number']
    run = {'id': tmp['id']}
    for a in _attrs:
        run[a] = cnv[a]
    #print(run)
    return run

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

# def getCol_Stats(runs, attrib):
#     """Returns the column values over all of the samples"""
#     ls = []
#     for r in runs:
#         ls.append((r['tumor']['id'], r['tumor'][attrib]))
#         ls.append((r['normal']['id'], r['normal'][attrib]))
#     vals = [x[1] for x in ls]
#     mean = sum(vals)/len(vals)
#     std = np.std(vals)
#     return (ls, mean, std)

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = map(lambda x: prettyprint(x), _attrs)
    hdr = ['Run'] + list(hdr)
    out.write("%s\n" % ",".join(hdr))
    for r in runs:
        tmp = [r['id']]
        for a in _attrs:
            tmp.append(str(r[a]))
        out.write("%s\n" % ",".join(tmp))
    out.close()

if __name__ == '__main__':
    main()
