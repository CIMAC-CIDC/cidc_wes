#!/usr/bin/env python
"""Len Taing 2020 (TGBTG)"""

import os
import sys
import subprocess
import json
import base64

from optparse import OptionParser

def parseFile(in_file):
    """parse the relevant file--edit this!"""
    f = open(in_file)
    hdr = f.readline().strip().split()
    l = f.readline().strip().split()
    tmp = dict(zip(hdr, l))

    ret = {'purity': tmp['purity'],
           'plodiy': tmp['ploidy'],
           'dipLogR': tmp['dipLogR']}
    #print(ret)
    f.close()
    return ret

def main():
    usage = "USAGE: %prog -r run_name -f file -o output_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-f", "--file", help="purity file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    js_out = {'id': options.run}
    purity = parseFile(options.file)
    js_out['copy_number'] = purity

    out = open(options.output, 'w')
    out.write(json.dumps(js_out))
    out.close()
    
if __name__=='__main__':
    main()
