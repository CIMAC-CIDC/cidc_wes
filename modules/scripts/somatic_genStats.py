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
    optparser.add_option("-o", "--out", help="output mutation variant type count .csv file")
    optparser.add_option("-n", "--nonsyn", help="output nonsynonymous mutation count .csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.maf or not options.out or not options.nonsyn:
        optparser.print_help()
        sys.exit(-1)

    #READ in template file:
    cts = {}
    non_cts = {}
    for maf_f in options.maf:
        f = open(maf_f)
        #HACK: get run name from file
        run_name = maf_f.split("/")[-1].split("_")[0]
        #print(run_name)
        run_ct = {"SNP": 0, "INS": 0, "DEL": 0}
        non_syn_ct = {"SNP": 0, "INS": 0, "DEL": 0}
        for l in f:
            if l.startswith("#"):
                continue
            tmp = l.strip().split("\t")
            #variant classification is 9th col, e.g. Missense_Mutation, Nonsense_Mutation
            #variant type is 10th col, e.g. SNP/INS/DEL
            classification = tmp[8]
            label = tmp[9]
            if label in run_ct:
                #print(label)
                run_ct[label] += 1
            #HANDLE non-synomous i.e. Missense_Mutation or Nonsense_Mutation
            if classification in ['Missense_Mutation', 'Nonsense_Mutation']:
                if label in non_syn_ct:
                    non_syn_ct[label] += 1
        f.close()
        run_ct['TOTAL']= sum([run_ct['SNP'], run_ct['INS'], run_ct['DEL']])
        non_syn_ct['TOTAL']= sum([non_syn_ct['SNP'], non_syn_ct['INS'], non_syn_ct['DEL']])
        
        cts[run_name] = run_ct
        non_cts[run_name] = non_syn_ct

    out = open(options.out, "w")
    non_syn_out = open(options.nonsyn, "w")
    out.write("%s\n" % ",".join(["Run","TOTAL", "SNP","INSERT","DELETE"]))
    non_syn_out.write("%s\n" % ",".join(["Run","TOTAL", "SNP","INSERT","DELETE"]))
    for r in sorted(cts.keys()):
        tmp = cts[r]
        tmp_non = non_cts[r]
        out.write("%s\n" % ",".join([r, str(tmp['TOTAL']), str(tmp["SNP"]),
                                     str(tmp["INS"]), str(tmp["DEL"])]))
        non_syn_out.write("%s\n" % ",".join([r, str(tmp_non['TOTAL']),
                                             str(tmp_non["SNP"]),
                                             str(tmp_non["INS"]),
                                             str(tmp_non["DEL"])]))
    out.close()
    non_syn_out.close()

if __name__=='__main__':
    main()
