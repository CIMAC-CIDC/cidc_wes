#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import csv
import numpy as np
from optparse import OptionParser

#Attribs to read in
#LEAVE off the plot for now
#_attrs = ['mutation_summary', 'transition_matrix', 'tmb', 'functional_summary']

def processJson(json_fpath, annotations):
    """Given a json file, returns 1. the entire json record and  
    2. a dictionary of that file with the values needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make runs
    run = {'id': tmp['id'], 'tumor': tmp['tumor']['id'], 'normal':tmp['normal']['id']}
    #Add annotations
    ann = annotations.get(tmp['id'])
    run['annotations'] = {}
    if ann:
        for a in ann:
            run['annotations'][a] = ann[a]

    #print(run)
    return (tmp, run)

def parseMeta(metasheet):
    """Parse the metasheet.cohort.csv file, skipping comments
    https://stackoverflow.com/questions/14158868/python-skip-comment-lines-marked-with-in-csv-dictreader
    """
    annotations = {}
    metasheet_f = open(metasheet)
    reader = csv.DictReader(filter(lambda row: row[0]!='#', metasheet_f))
    for row in reader:
        run = row["Run"]
        del row["Run"] #The run col is not an annotation
        annotations[run] = {} #init the dict
        for ann in row.keys():
            annotations[run][ann] = row[ann]
    metasheet_f.close()
    #print(annotations)
    return annotations

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-m", "--metasheet", help="metasheet.cohort.csv")
    optparser.add_option("-r", "--runs_output", help="json output - runs in json form")
    optparser.add_option("-d", "--data_output", help="json output - sample json data")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.metasheet or not options.runs_output or not options.data_output:
        optparser.print_help()
        sys.exit(-1)

    #PARSE metasheet
    annotations = parseMeta(options.metasheet)

    #read in json data
    #NOTE: using zip(*) to unpack the list of tuples
    data, runs = zip(*[processJson(f, annotations) for f in options.files])
    #print(runs)
    #print(data)

    #write runs output
    json_out = open(options.runs_output, "w")
    json_out.write(json.dumps(runs, indent=4))
    json_out.close()

    #write data output
    json_out = open(options.data_output, "w")
    json_out.write(json.dumps(data, indent=4))
    json_out.close()

    #DEPRECATED!
    # #Samples view:
    # samples = []
    # for r in runs:
    #     tmr = {'id': r['tumor'], 'tissue':'Tumor'}
    #     nrm = {'id': r['normal'], 'tissue':'Normal'}
    #     for k in r.keys():
    #         if k != 'id' and k != 'tumor' and k != 'normal':
    #             tmr[k] = r[k]
    #             nrm[k] = r[k]
    #     samples.append(tmr)
    #     samples.append(nrm)
    # json_out = open(options.samples_output, "w")
    # json_out.write(json.dumps(samples, indent=4))
    # json_out.close()

if __name__ == '__main__':
    main()
