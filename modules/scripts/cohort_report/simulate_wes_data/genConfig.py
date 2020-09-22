#!/usr/bin/env python3
"""Script to generate wes json data"""

import os
import sys

import random
import json
#import ruamel.yaml
import yaml

from optparse import OptionParser

def main():
    usage = "USAGE: %prog -d [dir] -o [output json file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-d", "--dir", help="directory to parse")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)


    if not options.dir:
        optparser.print_help()
        sys.exit(-1)

    tmp = {'runs':{}}
    files = os.listdir(options.dir)
    if options.output:
        out = open(options.output,"w")
    else:
        out = open("config.cohort.yaml", "w")
        
    for f in files:
        if f.endswith(".json"):
            fname = ".".join(f.split(".")[:-1])
            tmp['runs'][fname] = [os.path.abspath(os.path.join(options.dir, f))]
    #ruamel.yaml.round_trip_dump(wes_config, out)
    yaml.dump(tmp, out)
    out.close()
    
if __name__ == '__main__':
    main()
