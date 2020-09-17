#!/usr/bin/env python3
"""
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [filtered.tsv] -o [output file.txt]"
    optparser = OptionParser(usage=usage)                    
    optparser.add_option("-f", "--file", help="neoantigen filter.tsv file")
    optparser.add_option("-o", "--out", help="output file")

    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.file or not options.out:
        optparser.print_help()
        sys.exit(-1)
        
    f = open(options.file)
    #Fields we want to pick out from the pvacseq filtered output file
    fields = ['Gene Name', 'Mutation', 'Protein Position', 'HGVSc', 'HGVSp',
              'HLA Allele', 'MT Epitope Seq', 'MT IC50', 'WT IC50',
              'Fold Change', 'Tumor DNA VAF', 'Score']
    hdr = f.readline().strip().split("\t")
    results = []
    for l in f:
        tmp = dict(zip(hdr, l.strip().split("\t")))
        ls = [tmp[fld] for fld in fields]
        #print(ls)
        results.append(ls)
    f.close()

    out = open(options.out, 'w')
    #print out the header
    out.write("%s\n" % "\t".join(fields))
    for row in results:
        out.write("%s\n" % "\t".join(row))
    out.close()
    
if __name__=='__main__':
    main()
