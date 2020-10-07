#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

def removePath(image_f, report_path):
    return "img:%s" % image_f.replace(report_path, "")

def main():
    usage = "USAGE: %prog -f [lollipop plot png file path] -f [lollipop plot png file path] -r [cohort report path] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="lollipop png file path")
    optparser.add_option("-r", "--path", help="report path", default="")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #Ensure that the report path has a trailing /
    if options.path and not options.path.endswith("/"):
        options.path = options.path + "/"

    #MAIN output.tsv
    _group = 3
    out = open(options.output, "w")
    hdr = ["&nbsp;"] * (_group+1) #an invisible hdr
    out.write("%s\n" % ",".join(hdr))
    
    #group images into by 3
    nrows = math.ceil(len(options.files)/float(_group))
    for i in range(nrows):
        start= i*_group
        end = start + _group
        if end > len(options.files):
            end = len(options.files)
        row = ["&nbsp;"] #make first col invisible
        tmp = [removePath(f, options.path) for f in options.files[start:end]]
        row.extend(tmp)
        out.write("%s\n" % ",".join(row))
    out.close()

if __name__ == '__main__':
    main()
