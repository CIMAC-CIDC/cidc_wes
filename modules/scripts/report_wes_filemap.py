#!/usr/bin/env python
"""
reads in the input files/yaml files and generates a yaml file for the entire 
wes run
"""

import os
import sys
import yaml
from string import Template
from optparse import OptionParser
            
def main():
    usage = "USAGE: %prog -v [wes_version_txt] -m [metrics all summary file] -y [list of file paths]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-v", "--wes_version", help="wes_version.txt path")
    optparser.add_option("-m", "--metrics_summary", help="metrics all summary file path")
    optparser.add_option("-y", "--files", action="append", help="list of output.yaml files")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.wes_version or not options.metrics_summary or not options.files:
        optparser.print_help()
        sys.exit(-1)

    samples = {}
    runs = {}

    for f in options.files:
        ffile = open(f)
        data = yaml.load(ffile)
        if 'samples' in data:
            for s in data['samples']:
                if not s in samples:
                    #init a new entry
                    samples[s]={}
                
                for k,v in data['samples'][s].items():
                    samples[s][k] = v
        elif 'runs' in data:
            for r in data['runs']:
                if not r in runs:
                    #init a new entry
                    runs[r]={}
                
                for k,v in data['runs'][r].items():
                    runs[r][k] = v

        ffile.close()
    #print(samples)
    #print(runs)

    out = {'wes_version': options.wes_version,
           'metrics_all_sample_summaries': options.metrics_summary,
           'samples': samples,
           'runs': runs}
    print(yaml.dump(out))

if __name__=='__main__':
    main()
