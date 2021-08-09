#!/usr/bin/env python
"""Len Taing 2021 (TGBTG)
Given a version of the sample summary json 
AND the {run}.filter.tsv file this script will: 

Add a section 'neoantigen' to the json file.  The json section will have the 
format: {'neoantigen_filter_file': <base64 encoding of entire filter.tsv file>}

Output: saves the resulting json file to the specified output fname path

NOTE: Many of these fns are taken directly from json_neoantigen.py 
"""

import os
import io
import sys
import subprocess
import json
import base64
from optparse import OptionParser

#KEY component: reusing the key functions directly from mutProfile.py
#Assume that THIS script resides in cidc_wes/shell_scripts, and we want to
#inject cidc_wes/cidc-vs

#Try to find the first occurance of 'cidc_wes' drops everything after
_index=(os.path.abspath(__file__).split("/").index('cidc_wes'))
_cidc_wes_path=os.path.sep.join(os.path.abspath(__file__).split("/")[:_index+1])
sys.path.append(os.path.join(_cidc_wes_path, 'modules', 'scripts'))
#print(sys.path)

#import somatic_getTopOncoGenes

def main():
    usage = "USAGE: %prog -j sample summary file (json) -r [ABSOLUTE PATH to .fasta reference] -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-j", "--json_file", help="sample summary json file", default=None)
    optparser.add_option("-f", "--filter_file", help="run_filter.tsv file from pvacseq")
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.json_file or not options.filter_file or  not options.output:
        optparser.print_help()
        sys.exit(-1)

    #READ in the filter_file and convert to base64 string
    f = open(options.filter_file)
    s = f.read()
    s_byte = s.encode('utf-8')
    #NOTE: .decode is needed to convert the bytes back to a string
    s_b64 = base64.b64encode(s_byte).decode('utf-8')
    f.close()

    #READ in the sample json file
    f = open(options.json_file)
    sample_json = json.load(f)
    f.close()

    neoantigen = {'neoantigen_filter_file': s_b64}
    
    #integrate this back into sample_json
    sample_json['neoantigen'] = neoantigen

    #DUMP the output
    out = open(options.output, 'w')
    json.dump(sample_json, out)
    out.close()
    
if __name__=='__main__':
    main()
