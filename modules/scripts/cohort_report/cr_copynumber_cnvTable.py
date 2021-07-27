#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

from cr_utils import prettyprint
#Attribs to read in
#LEAVE off the plot for now
#_attrs = ['clonality', 'purity', 'ploidy', 'dipLogR']#, 'cnv_plot_file', 'cnv_file']

def processJson(json_fpath, attrs):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    run = None
    #make samples
    if 'copy_number' in tmp:
        cnv = tmp['copy_number']
        run = {'id': tmp['id']}
        for a in attrs:
            run[a] = cnv[a]
        #print(run)
    return run

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-a", "--attributes", action="append", help="attributes to select/output")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.attributes or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #try to get attributes list
    _attrs = options.attributes
    
    #read in json data
    runs = [processJson(f, _attrs) for f in options.files]
    #NEED to remove samples without 'copy_number' section, ie tumor-only sample
    runs = filter(lambda x: x, runs)

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
