#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
import gzip
import subprocess
import base64
from optparse import OptionParser

#Attribs to read in
#LEAVE off the plot for now
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
    
def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output gz file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    out = gzip.open(options.output, 'wb')
    for maf in [r['filtered_maf_file'] for r in runs]:
        #NOTE: we don't need to add .decode('utf-8') to end of base64 call b/c
        #GZIP takes in bytes, not strings
        out.write(base64.b64decode(maf))
    out.close()

if __name__ == '__main__':
    main()
