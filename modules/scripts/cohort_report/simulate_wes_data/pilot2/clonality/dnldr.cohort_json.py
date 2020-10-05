#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import yaml
from optparse import OptionParser

def parseYaml(yaml_file):
    """Parses a yaml file and returns a dictionary"""
    ret = {}
    with open(yaml_file, 'r') as stream:
        try:
            ret = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return ret

def main():
    usage = "USAGE: %prog -c config.cohort_json.yaml -f [input file of google bucket paths]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--config", help="config yaml file that defines runs and samples", default=None)
    optparser.add_option("-f", "--file", help="input file of google bucket paths", default=None)
    optparser.add_option("-o", "--out_dir", help="directory to download the files", default=".")
    optparser.add_option("-r", "--runs", help="Runs level files", action="store_true")
    optparser.add_option("-s", "--subdirs", help="make subdirs in the out_dir for each run/sample", action="store_true")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.config or not options.file: # or not options.ext:
        optparser.print_help()
        sys.exit(-1)
    
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    config = parseYaml(options.config)
    runs = list(config.keys())
    samples = []
    for r in runs:
        samples.append(config[r]['tumor'])
        samples.append(config[r]['normal'])

    added_list= [] #To ensure that we don't duplicate, check against this list
    dnld_list = []
    f = open(options.file)
    check_list = runs if options.runs else samples #Check if this is run/sample level
    for l in f:
        path = l.strip()
        #KEY in on the sub-dir
        #NOTE: we are not keying in on the file name b/c we assume wes_results
        #was built using a grep cmd for that
        subdir = path.split("/")[-3]
        #print(subdir)
        if subdir in check_list and subdir not in added_list:
            dnld_list.append(path)
            added_list.append(subdir)
            out_path = options.out_dir
            if options.subdirs:
                out_path = os.path.join(options.out_dir, subdir)
            #print(out_path)
            #create the out_path if it's not existent
            if not os.path.isdir(out_path):
                os.mkdir(out_path)
            #Download
            cmd = "gsutil -m cp %s %s" % (path, out_path)
            cmd = cmd.split(" ")
            output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    #print(dnld_list)
    f.close()
    
if __name__=='__main__':
    main()
