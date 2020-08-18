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
_attrs = ['mutation_summary']

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

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-o", "--output", help="output tsv file")
    optparser.add_option("-j", "--json_output", help="json output - mutation summary table in json form")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.output or not options.json_output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #build mutation summaries dict
    mutation_summaries = {}

    ###########################################################################
    #KEY: need to make this tsv file b/c otherwise the arguments sent to
    #     the jscript fn will be read as cols
    ###########################################################################
    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = ['Run', 'Total Mutations']
    out.write("%s\n" % "\t".join(hdr))
    for r in runs:
        tmp = [r['id']]

        #handle mutation summary - Show total
        total = str(r['mutation_summary']['total'])
        #tmp.append("""[%s](#){: onclick=\"wes_modal_helper(somatic_summary, %s, table);\" }""" % (total, r['id']))
        tmp.append(total)
        #tmp.append("<a href='http://www.yahoo.com'>ref1</a>")
        out.write("%s\n" % "\t".join(tmp))

        mutation_summaries[r['id']] = r['mutation_summary']
    out.close()

    json_out = open(options.json_output, "w")
    json_out.write(json.dumps(mutation_summaries))
    json_out.close()

if __name__ == '__main__':
    main()
