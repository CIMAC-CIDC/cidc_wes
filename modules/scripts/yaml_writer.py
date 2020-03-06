#!/usr/bin/env python
"""
simple yaml writer--takes input args, key and value lists and dumps them to 
the std out as valid yaml
"""

import os
import sys
import yaml
from string import Template
from optparse import OptionParser
            
def main():
    usage = "USAGE: %prog -t [title either samples or runs] -k [list of key names] -f [list of file paths]\nNOTE: key list and file list lengths must match!"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-t", "--title", help="either samples or runs")
    optparser.add_option("-n", "--name", help="name of sample/run")
    optparser.add_option("-k", "--keys", action="append", help="list of key names")
    optparser.add_option("-f", "--files", action="append", help="list of file paths")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.title or not options.name or not options.keys or not options.files:
        optparser.print_help()
        sys.exit(-1)

    tmp = {options.title: {options.name: dict([(k,options.files[i]) for (i, k) in enumerate(options.keys)])}}
    print(yaml.dump(tmp))

if __name__=='__main__':
    main()
