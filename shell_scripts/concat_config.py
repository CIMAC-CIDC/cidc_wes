#!/usr/bin/env python
"""Len Taing 2021 (TGBTG) Jacob Geisberg 2022
Given a CSV file of tumor,normal will generate a yaml file to input into
concat.snakefile
"""

import os
import sys
import subprocess
import yaml
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f csv file of t-n pairings"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="csv file of t-n pairings", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)
    optparser.add_option("-s", "--source", help="path to bucket where fastqs are stored", default=None)
    optparser.add_option("-d", "--destination", help="path for merged fastq files to be sent", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file or not options.output or not options.source or not options.destination:
        optparser.print_help()
        sys.exit(-1)

    path = options.source
    destination = options.destination
    tmp = {'path': path, 'destination':destination, 'samples': {}}
    with open(options.file) as f:
        for l in f:
            (tmr, nrm) = l.strip().split(",")
            #Add tumor pairs
            tmr_r1 = ["%s/%s/r1_L1.fastq.gz" % (path,tmr),
                      "%s/%s/r1_L2.fastq.gz" % (path,tmr)]
            tmr_r2 = ["%s/%s/r2_L1.fastq.gz" % (path,tmr),
                      "%s/%s/r2_L2.fastq.gz" % (path,tmr)]
            tmp['samples']['%s_R1' % tmr] = tmr_r1
            tmp['samples']['%s_R2' % tmr] = tmr_r2
            #Add normal pairs
            nrm_r1 = ["%s/%s/r1_L1.fastq.gz" % (path,nrm),
                      "%s/%s/r1_L2.fastq.gz" % (path,nrm)]
            nrm_r2 = ["%s/%s/r2_L1.fastq.gz" % (path,nrm),
                      "%s/%s/r2_L2.fastq.gz" % (path,nrm)]
            tmp['samples']['%s_R1' % nrm] = nrm_r1
            tmp['samples']['%s_R2' % nrm] = nrm_r2
    
    yaml_out = open(options.output, "w")
    yaml.dump(tmp, yaml_out) 
    yaml_out.close()
    
if __name__=='__main__':
    main()
