#!/usr/bin/env python
"""Jacob Geisberg 2022"""
import pandas as pd
from optparse import OptionParser
import sys
import subprocess #not sure this is needed

def main():
    usage = "USAGE: %prog -r run_name -f file -o output_file"
    optparser = OptionParser(usage=usage)
    #optparser.add_option("-r", "--run", help="run name", default=None)
    optparser.add_option("-f", "--file", help="tcellextrect file", default=None)
    optparser.add_option("-o", "--output", help="output file", default=None)

    (options, args) = optparser.parse_args(sys.argv)

    if not options.file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    df = pd.read_csv(options.file)
    df = df[['sample','TCRA.tcell.fraction', 'qcFit']]
    df = df.rename(columns={"sample": "Run",
                            "TCRA.tcell.fraction": "Tcell_fraction",
                            "qcFit": "q_value"})
    df = df.round(4)
    #df['Run'][0] = options.run
    df.to_csv(options.output, index=False)


if __name__=='__main__':
    main()


          
