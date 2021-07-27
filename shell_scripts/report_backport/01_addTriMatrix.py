#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)
Given a version of the sample summary report that includes the full MAF file
from somatic variant calling, this script will calculate the trinucleotide 
matrix and inject it into the json file

Output: saves the resulting json file to the specified output fname path

NOTE: Many of these fns are taken directly from mutProfile.py written by 
Jingxin Fu
NOTE: mutProfile.py requires the ABSOLUTE PATH to the references!!!
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
sys.path.append(os.path.join(_cidc_wes_path, 'cidc-vs'))
#print(sys.path)
import mutProfile

def main():
    usage = "USAGE: %prog -j sample summary file (json) -r [ABSOLUTE PATH to .fasta reference] -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-j", "--json_file", help="sample summary json file", default=None)
    optparser.add_option("-r", "--ref", help="ABSOLUTE path to fasta reference file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.json_file or not options.ref or  not options.output:
        optparser.print_help()
        sys.exit(-1)

    #READ in the sample json file
    f = open(options.json_file)
    sample_json = json.load(f)
    f.close()

    # decode somatic.filtered_maf_file
    maf_bytes = base64.b64decode(sample_json['somatic']['filtered_maf_file'])
    maf = maf_bytes.decode('utf-8')

    #Need to convert the string to a STREAM for pandas read-in
    maf_io = io.StringIO(maf)
    tri_mtrx_df = mutProfile.genTrinucleotideMtrx(maf_io, options.ref)
    tri_mtrx_dict = mutProfile.convertTriMatrix(tri_mtrx_df)

    #integrate this back into sample_json
    sample_json['somatic']['tri_matrix'] = tri_mtrx_dict

    #DUMP the output
    out = open(options.output, 'w')
    json.dump(sample_json, out)
    out.close()
    
if __name__=='__main__':
    main()
