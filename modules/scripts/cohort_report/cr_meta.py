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
#_attrs = ['mutation_summary', 'transition_matrix', 'tmb', 'functional_summary']

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    foo = tmp['meta']
    #print(tmp)
    run = {'id': tmp['id'], 'tumor': tmp['tumor']['id'], 'normal':tmp['normal']['id']}
    for a in foo.keys():
        run[a] = foo[a]
    #print(run)
    #sys.exit()
    return run
    
def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-r", "--runs_output", help="json output - runs in json form")
    optparser.add_option("-s", "--samples_output", help="json output - samples in json form")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.runs_output or not options.samples_output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #write runs output
    json_out = open(options.runs_output, "w")
    json_out.write(json.dumps(runs))
    json_out.close()

    #Samples view:
    samples = []
    for r in runs:
        tmr = {'id': r['tumor'], 'tissue':'Tumor'}
        nrm = {'id': r['normal'], 'tissue':'Normal'}
        for k in r.keys():
            if k != 'id' and k != 'tumor' and k != 'normal':
                tmr[k] = r[k]
                nrm[k] = r[k]
        samples.append(tmr)
        samples.append(nrm)
    json_out = open(options.samples_output, "w")
    json_out.write(json.dumps(samples))
    json_out.close()

if __name__ == '__main__':
    main()
