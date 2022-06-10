#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
import argparse


def integrate(dict_a, dict_omni):
    """Given an omnibus dict , we want to integrate elements of dict_a into
    it without over-writing existing elements.
    Recursively calls when a sub-elemtn is a dictionary
    """
    for k in dict_a.keys():
        #check is k exists
        if not k in dict_omni: #assign element
            dict_omni[k] = dict_a[k]
        else:
            #Conflict--if it's a dictionary element recur
            if isinstance(dict_omni[k],dict):
                integrate(dict_a[k], dict_omni[k])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--run", required=True, help="run name", default=None)
    parser.add_argument("-m", "--mapping", required=True, help="mapping json file", default=None)
    parser.add_argument("-c", "--coverage", required=True, help="coverage json file", default=None)
    parser.add_argument("-g", "--gc", required=True, help="gc content json file", default=None)
    parser.add_argument("-i", "--insert_size", required=True, help="insert size json file", default=None)
    parser.add_argument("-q", "--mean_quality", required=True, help="mean quality json file", default=None)
    parser.add_argument("-j", "--hla", required=True, help="hla json file", default=None)
    parser.add_argument("-p", "--purity", help="purity json file", default=None)
    parser.add_argument("-s", "--somatic", required=True, help="somatic filtered maf json file", default=None)
    parser.add_argument("-l", "--clonality", help="clonality json file", default=None)
    parser.add_argument("-n", "--neoantigen", help="neoantigen json file", default=None)
    parser.add_argument("-t", "--tcellextrect", help="tcellextrect json file", default=None)
    parser.add_argument("-e", "--msisensor2", help="msisensor2 json file", default=None)
    parser.add_argument("-v", "--copynumber", help="copynumber json file", default=None)
    parser.add_argument("-o", "--output", required=True, help="output file", default=None)

    # parse the argumnets and convert namespace object to dictionary
    args = parser.parse_args()
    args = vars(args)

    # extract the run and output arguments from the dictionary to process them separately
    run = args.pop('run') 
    output = args.pop('output')

    # get list of json files from args and process them one by one
    files = [args[value] for value in args if value is not None]
    js_out = {'id': run, 'meta': {}}
    for json_f in files:

        if json_f:
            #read in json
            f = open(json_f)
            tmp = json.load(f)
            f.close()

            #integrate tmp dict into js_out
            integrate(tmp, js_out)


    json_out = open(output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()



if __name__=='__main__':
    main()
