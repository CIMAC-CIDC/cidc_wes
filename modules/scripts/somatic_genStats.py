#!/usr/bin/env python
"""
generate somatic mutation statistics. output .csv file
"""

import os
import sys
from string import Template
from optparse import OptionParser
            
def main():
    usage = "USAGE: %prog -m [filter.maf files list] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--maf", action="append", help="filter.maf files")
    optparser.add_option("-o", "--out", help="output .csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.maf or not options.out:
        optparser.print_help()
        sys.exit(-1)

    #READ in template file:
    cts = {}
    for maf_f in options.maf:
        f = open(maf_f)
        #HACK: get run name from file
        run_name = maf_f.split("/")[-1].split("_")[0]
        #print(run_name)
        run_ct = {"SNP": 0, "INS": 0, "DEL": 0}
        for l in f:
            if l.startswith("#"):
                continue
            tmp = l.strip().split("\t")
            #Classification is 10th col
            label = tmp[9]
            if label in run_ct:
                #print(label)
                run_ct[label] += 1
        f.close()
        run_ct['TOTAL']= sum([run_ct['SNP'], run_ct['INS'], run_ct['DEL']])
        cts[run_name] = run_ct

    out = open(options.out, "w")
    out.write("%s\n" % ",".join(["Run","TOTAL", "SNP","INSERT","DELETE"]))
    for r in sorted(cts.keys()):
        tmp = cts[r]
        out.write("%s\n" % ",".join([r, str(tmp['TOTAL']), str(tmp["SNP"]),
                                     str(tmp["INS"]), str(tmp["DEL"])]))
    out.close()

if __name__=='__main__':
    main()
