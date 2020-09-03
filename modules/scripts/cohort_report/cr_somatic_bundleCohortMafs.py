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
    optparser.add_option("-o", "--output", help="output gz file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    out = gzip.open(options.output, 'wb')
    for maf in [r['filtered_maf_file'] for r in runs]:
        f = open(maf)
        out.write(bytes(f.read(), 'utf-8'))
        f.close()
    out.close()
    # mafs = [r['filtered_maf_file'] for r in runs]
    # cmd = "cat %s | gzip > %s" % (" ".join(mafs), options.output)
    # #cmd = "gzip -c %s > %s" % (" ".join(mafs), options.output)
    # #subprocess.run(cmd.split(" "))
    # output, err = subprocess.Popen(cmd.split(" "), #stdout=subprocess.PIPE,
    #                                stderr=subprocess.PIPE).communicate()

if __name__ == '__main__':
    main()
