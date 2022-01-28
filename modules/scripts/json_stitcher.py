#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser


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
    usage = "USAGE: %prog -r run_name ...."
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-m", "--mapping", help="mapping json file", default=None)
    optparser.add_option("-c", "--coverage", help="coverage json file", default=None)
    optparser.add_option("-g", "--gc", help="gc content json file", default=None)
    optparser.add_option("-i", "--insert_size", help="insert size json file", default=None)
    optparser.add_option("-q", "--mean_quality", help="mean quality json file", default=None)
    optparser.add_option("-j", "--hla", help="hla json file", default=None)
    optparser.add_option("-p", "--purity", help="purity json file", default=None)
    optparser.add_option("-s", "--somatic", help="somatic filtered maf json file", default=None)
    optparser.add_option("-l", "--clonality", help="clonality json file", default=None)
    optparser.add_option("-n", "--neoantigen", help="neoantigen json file", default=None)
    optparser.add_option("-t", "--tcellextrect", help="tcellextrect json file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.mapping or not options.coverage or not options.gc or not options.insert_size or not options.mean_quality or not options.hla or not options.somatic or not options.neoantigen or not options.tcellextrect or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #missing meta information
    js_out = {'id': options.run, 'meta': {}}
    for json_f in [options.mapping, options.coverage, options.gc,options.insert_size,
                   options.mean_quality, options.hla,options.purity, options.somatic,
                   options.clonality, options.neoantigen, options.tcellextrect]:

        if json_f:
            #read in json
            f = open(json_f)
            tmp = json.load(f)
            f.close()

            #integrate tmp dict into js_out
            integrate(tmp, js_out)
        
    
    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
