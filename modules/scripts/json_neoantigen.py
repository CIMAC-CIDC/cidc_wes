#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
import base64
from optparse import OptionParser


def main():
    usage = "USAGE: %prog -f file -o output json file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-f", "--in_file", help="input maf file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.in_file or  not options.output:
        optparser.print_help()
        sys.exit(-1)


    f = open(options.in_file)
    s = f.read()
    s_byte = s.encode('utf-8')
    #NOTE: .decode is needed to convert the bytes back to a string
    s_b64 = base64.b64encode(s_byte).decode('utf-8')
    f.close()

    js_out = {'id': options.run, 'neoantigen': {'neoantigen_filter_file':"%s" % s_b64}}

    out = open(options.output, 'w')
    out.write(json.dumps(js_out))
    out.close()
    
if __name__=='__main__':
    main()
