#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

#Attribs to read in
#LEAVE off the plot for now
_attrs = ['mutation_summary', 'transition_matrix', 'tmb', 'functional_summary']

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    somatic = tmp['somatic']
    #print(somatic)
    run = {'id': tmp['id']}
    for a in _attrs:
        run[a] = somatic[a]
    #print(run)
    #sys.exit()
    return run

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def calcTiTv(trans_mat):
    """Calculates the Transition/Transversion (Ti/Tv) ratio"""
    _transitions = ["AG", "GA", "CT", "TC"]
    total = 0
    ti_ct = 0
    for (k, row) in trans_mat.items():
        for (base, n) in row.items():
            total += n
            tmp = k + base
            if tmp in _transitions:
                ti_ct += n
    #print(total, ti_ct)
    ratio = "%.3f" % (float(ti_ct)/float(total-ti_ct))
    return ratio
    
def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output tsv file")
    optparser.add_option("-j", "--json_output1", help="json output - mutation summary table in json form")
    optparser.add_option("-k", "--json_output2", help="json output - transition matrix table in json form")
    optparser.add_option("-l", "--json_output3", help="json output - tmb table in json form")
    optparser.add_option("-m", "--json_output4", help="json output - functional summary table in json form")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output or not options.json_output1 or not options.json_output2 or not options.json_output3 or not options.json_output4:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #build mutation summaries dict
    mutation_summaries = {}
    #transition matrices dict
    trans_mat = {}
    tmb = {}
    functional_summaries = {}
    
    ###########################################################################
    #KEY: need to make this tsv file b/c otherwise the arguments sent to
    #     the jscript fn will be read as cols
    ###########################################################################
    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = ['Run', 'Total Mutations', 'TiTv', 'TMB', 'Nonsyn Mutations']
    out.write("%s\n" % "\t".join(hdr))
    for r in runs:
        tmp = [r['id']]

        #handle mutation summary - Show total
        total = str(r['mutation_summary']['total'])
        #tmp.append("""[%s](#){: onclick=\"wes_modal_helper(somatic_summary, %s, table);\" }""" % (total, r['id']))
        tmp.append(total)

        #Transition matrix
        ti_tv = calcTiTv(r['transition_matrix'])
        tmp.append(ti_tv)

        #TMB: take Tumor
        tmp.append(str(r['tmb']['tumor']))

        #Nonsyn. mutations
        missense = r['functional_summary']['missense']['total']
        nonsense = r['functional_summary']['nonsense']['total']
        tmp.append(str(missense + nonsense))
        
        #tmp.append("<a href='http://www.yahoo.com'>ref1</a>")
        out.write("%s\n" % "\t".join(tmp))

        mutation_summaries[r['id']] = r['mutation_summary']
        trans_mat[r['id']] = r['transition_matrix']
        tmb[r['id']] = r['tmb']
        functional_summaries[r['id']] = r['functional_summary']
        
    out.close()

    json_out = open(options.json_output1, "w")
    json_out.write(json.dumps(mutation_summaries))
    json_out.close()

    json_out = open(options.json_output2, "w")
    json_out.write(json.dumps(trans_mat))
    json_out.close()

    json_out = open(options.json_output3, "w")
    json_out.write(json.dumps(tmb))
    json_out.close()

    json_out = open(options.json_output4, "w")
    json_out.write(json.dumps(functional_summaries))
    json_out.close()

if __name__ == '__main__':
    main()
