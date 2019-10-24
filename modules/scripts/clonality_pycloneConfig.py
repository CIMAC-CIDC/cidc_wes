#!/usr/bin/env python
"""
generate pyclone config bassed off of a template
"""

import os
import sys
from string import Template
from optparse import OptionParser
            
def main():
    usage = "USAGE: %prog -t [pyclone config template] -m [mutation file path] -o [output dir]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-t", "--template", help="template file")
    optparser.add_option("-m", "--mutation_file", help="mutation file")
    optparser.add_option("-o", "--outdir", help="output dir")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.template or not options.mutation_file or not options.outdir:
        optparser.print_help()
        sys.exit(-1)

    #READ in template file:
    f = open(options.template)
    conf = Template(f.read())
    f.close()
    print(conf.substitute(mutation_file_full_path=options.mutation_file,
                          outdir_full_path=options.outdir))

if __name__=='__main__':
    main()
