#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)
Given a version of the sample summary json that includes the full MAF file
from somatic variant calling, AND a list of top onco driving genes, 
this script will: 

Add a section 'geneList' to the somatic section where each gene represents
a top onco gene that is mutated in the sample i.e. it's found in the maf
Output: saves the resulting json file to the specified output fname path

NOTE: Many of these fns are taken directly from somatic_getTopOncoGenes.py 
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
import somatic_getTopOncoGenes

def main():
    usage = "USAGE: %prog -j sample summary file (json) -r [ABSOLUTE PATH to .fasta reference] -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-j", "--json_file", help="sample summary json file", default=None)
    optparser.add_option("-l", "--oncoGeneList", help="list of top onco driving genes")
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.json_file or not options.oncoGeneList or  not options.output:
        optparser.print_help()
        sys.exit(-1)

    #READ in the cancerGeneList
    oncoGeneList = []
    f = open(options.oncoGeneList)
    for l in f:
        oncoGeneList.append(l.strip())
    f.close()

    #READ in the sample json file
    f = open(options.json_file)
    sample_json = json.load(f)
    f.close()

    # decode somatic.filtered_maf_file
    maf_bytes = base64.b64decode(sample_json['somatic']['filtered_maf_file'])
    maf = maf_bytes.decode('utf-8').strip().split("\n")

    #Calculate variant summaries
    geneList = somatic_getTopOncoGenes.getOncoGeneList(maf, oncoGeneList, True)
    
    #integrate this back into sample_json
    sample_json['somatic']['geneList'] = geneList

    #DUMP the output
    out = open(options.output, 'w')
    json.dump(sample_json, out)
    out.close()
    
if __name__=='__main__':
    main()
