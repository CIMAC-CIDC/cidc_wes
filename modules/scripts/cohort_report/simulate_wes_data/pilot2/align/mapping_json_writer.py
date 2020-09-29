#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
from optparse import OptionParser


def parseFile(map_file, dedup_file):
    """Reads in a _gc_metrics.txt file- 4th col"""
        #take 1st and 5th line from mapping file
    f = open(map_file)
    total = int(f.readline().strip().split(" ")[0])
    #skip next 4 lines
    for i in range(4):
        tmp = f.readline()
    mapped = int(f.readline().strip().split(" ")[0])
    f.close()

    #read in 5th line from dedup
    f = open(dedup_file)
    #skip first 4 lines
    for i in range(4):
        tmp = f.readline()    
    dedup = int(f.readline().strip().split(" ")[0])
    f.close()
    
    ret = {'total_reads': total, 'mapped_reads': mapped, 'dedup_reads': dedup}
    #print(ret)
    return ret

def main():
    usage = "USAGE: %prog -r run_name -t tumor_file -n normal_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-t", "--tumor_map", help="tumor sorted.bam mapping file", default=None)
    optparser.add_option("-s", "--tumor_dedup", help="tumor sorted.dedup.bam file", default=None)
    optparser.add_option("-n", "--normal_map", help="tumor sorted.bam mapping file", default=None)
    optparser.add_option("-m", "--normal_dedup", help="tumor sorted.dedup.bam file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.tumor_map or not options.tumor_dedup or not options.normal_map or not options.normal_dedup or not options.output:
        optparser.print_help()
        sys.exit(-1)

    js_out = {'id': options.run}
    tmr = parseFile(options.tumor_map, options.tumor_dedup)
    #GET tumor sample name
    fname = options.tumor_map.split("/")[-1]
    tmr_id = fname.split("_")[0]

    nrm = parseFile(options.normal_map, options.normal_dedup)
    #GET normal sample name
    fname = options.normal_map.split("/")[-1]
    nrm_id = fname.split("_")[0]

    #print(total, mapped, dedup)
    js_out['tumor'] = {'id': tmr_id, 'alignment': tmr}
    js_out['normal'] = {'id': nrm_id, 'alignment': nrm}
    
    json_out = open(options.output, "w")
    json_out.write(json.dumps(js_out))
    json_out.close()
    
if __name__=='__main__':
    main()
