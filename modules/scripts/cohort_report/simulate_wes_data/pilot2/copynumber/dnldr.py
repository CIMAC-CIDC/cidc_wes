#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [input file of google bucket paths]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="input file of google bucket paths", default=None)
    optparser.add_option("-o", "--out_dir", help="directory to download the files", default=".")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    f = open(options.file)
    for l in f:
        path = l.strip()
        cmd = "gsutil -m cp %s %s" % (path, options.out_dir)
        cmd = cmd.split(" ")
        output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    f.close()
    
if __name__=='__main__':
    main()
