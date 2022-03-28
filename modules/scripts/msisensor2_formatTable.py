#!/usr/bin/env python                                                                                            
"""Jacob Geisberg 2022"""
import pandas as pd
from optparse import OptionParser
import sys


def main():
    usage = "USAGE: %prog -r run_name -f file -o output_file"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="tcellextrect file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    df = pd.read_csv(options.file, sep='\t')
    df = df.round(4)
    df.columns = ['Total_Sites','Somatic_Sites','Percent_Somatic']
    df.to_csv(options.output, index=False)


if __name__=='__main__':
    main()



