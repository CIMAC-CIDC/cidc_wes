#!/usr/bin/env python
"""Len Taing (TGBTG) 2021

Given a Google bucket folder with this folder structure:
Example: gs://10021_results/WES/CNTZ36P8Q.01/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/report.tar.gz
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/align/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/copynumber/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/logs/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/metrics/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/msisensor2/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/neoantigen/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/optitype/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/report/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/rna/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/somatic/
   gs://10021_results/WES/CNTZ36P8Q.01/analysis/xhla/

Removes all of the folders under analysis EXCEPT for 'rna' and 'report.tar.gz'

The script uses an explicit list of sub-folders to delete
"""

import os
import sys
import subprocess
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -b [google bucket path, e.g. gs://10021_results/WES/CNTZ36P8Q.01/; NOTE: this folder is assumed to have an 'analysis' subfolder]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-b", "--bucket_path", help="google bucket path to the where the WES results are stored", default=None)
    optparser.add_option("-t", "--tumor_only", help="Is tumor-only sample", action="store_true", default=False)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.bucket_path:
        optparser.print_help()
        sys.exit(-1)
        
    if not options.bucket_path.startswith("gs://"):
        print("Does not appear to be a valid google bucket path %s" % options.bucket_path)
        sys.exit(-1)

    #Check if bucketpath ends with '/'
    bp = options.bucket_path
    if not options.bucket_path.endswith("/"):
        bp = bp + "/"

    analysis_path = bp + "analysis/"
    #BASE Sub folders to delete
    sub_folders = ['align', 'copynumber', 'metrics', 'msisensor2',
                   'neoantigen', 'optitype', 'report', 'somatic', 'xhla',
                   'hlahd', 'tcellextrect', 'cnvkit']
    
    if not options.tumor_only: #For Tumor-Normal samples add these
        sub_folders.extend(['clonality', 'corealignments', 'germline', 'purity'])

    for f in sub_folders:
        path = analysis_path + f
        cmd = "gsutil rm -r %s" % path
        print(cmd)
        output, err = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
if __name__=='__main__':
    main()
