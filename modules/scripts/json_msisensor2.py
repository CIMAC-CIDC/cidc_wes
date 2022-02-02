#!/usr/bin/env python                                                                                          
"""Jacob Geisberg 2022"""
import pandas as pd
import os
import sys
import subprocess
import json
from optparse import OptionParser

def parseFile(in_file):
    """parse the relevant file--edit this!"""

    f = pd.read_csv(in_file, sep='\t')
    f = f.round(4)
    ret = {'Total_Sites':int(f["Total_Number_of_Sites"][0]),
           'Somatic_Sites':int(f["Number_of_Somatic_Sites"][0]),
           'Percent_Somatic':int(f["%"][0])
          }


    return ret


def main():
    usage = "USAGE: %prog -r run_name -f file -o output_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-f", "--file", help="tcellextrect file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.run or not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    js_out = {'id': options.run}
    msisensor2 = parseFile(options.file)
    js_out['msisensor2'] = msisensor2

    out = open(options.output, 'w')
    out.write(json.dumps(js_out))
    out.close()

if __name__=='__main__':
    main()
