#!/usr/bin/env python
"""Len Taing 2021 (TGBTG)
Given a version of the sample summary json that includes the full MAF file
from somatic variant calling, this script will calculate the number of 
Total, SNP, Insert, DELETE, and TMB and inject it into the json file

Output: saves the resulting json file to the specified output fname path

NOTE: Many of these fns are taken directly from somatic_genStats.py 
NOTE: the targeted regions bed file is required to calculate TMB
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
import somatic_genStats

def main():
    usage = "USAGE: %prog -j sample summary file (json) -r [ABSOLUTE PATH to .fasta reference] -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-j", "--json_file", help="sample summary json file", default=None)
    optparser.add_option("-t", "--target", help="target region bed file")
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.json_file or not options.target or  not options.output:
        optparser.print_help()
        sys.exit(-1)

    #READ in the sample json file
    f = open(options.json_file)
    sample_json = json.load(f)
    f.close()

    # decode somatic.filtered_maf_file
    maf_bytes = base64.b64decode(sample_json['somatic']['filtered_maf_file'])
    maf = maf_bytes.decode('utf-8').strip().split("\n")

    #Calculate variant summaries
    (cts, annot_cts) = somatic_genStats.getVariantStats(maf, sample_json['id'], True)
    tmp = cts[sample_json['id']]
    
    #ADD tmb
    #Calculate the number of target_Mbs
    target_mb = round(somatic_genStats.countBedBases(options.target)/float(10**6), 1)
    TMB = round(tmp['TOTAL']/target_mb, 1)
    
    summary = {'total':tmp['TOTAL'], 'snp': tmp['SNP'], 'insert': tmp['INS'],
               'delete':tmp['DEL'], 'tmb': TMB}
    
    #integrate this back into sample_json
    sample_json['somatic']['summary'] = summary

    #DUMP the output
    out = open(options.output, 'w')
    json.dump(sample_json, out)
    out.close()
    
if __name__=='__main__':
    main()
