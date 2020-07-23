#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to generate wes run plots information"""

import os
import sys
import math
import json
import numpy as np
from optparse import OptionParser

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the 
    first letter of each word"""
    s = s.replace("_"," ")
    s = s.upper() if toUpper else s.title()
    return s

def getSampleInfo(d, isTumor):
    """Given a dictionary of json records, returns a record of tumor OR
    normal sample"""
    sample = d['tumor'] if isTumor else d['normal']
    alignment = sample['alignment']
    attrs = ['total_reads', 'mapped_reads', 'dedup_reads', 'mean_quality_score']
    tmp = {'id': sample['id']}
    for a in attrs:
        #millify at printout
        #tmp[a] = millify(alignment.get(a))
        tmp[a] = alignment.get(a)
    return tmp

def processJson(json_fpath):
    """Given a json file, returns a dictionary of that file with the values
    needed for this script"""
    #READ in json
    f = open(json_fpath)
    tmp = json.load(f)
    f.close()

    #make samples
    tumor = getSampleInfo(tmp, True)
    normal = getSampleInfo(tmp, False)
    run = {'id': tmp['id'], 'tumor': tumor, 'normal': normal}
    return run

def millify(n):
    """Given a large int n, returns a string representation of n in human-
    readable form
    ref: https://stackoverflow.com/questions/3154460/python-human-readable-large-numbers
    """
    millnames = ['',' K',' M',' B',' T']

    n = float(n)
    millidx = max(0,min(len(millnames)-1,
                    int(math.floor(0 if n == 0 else math.log10(abs(n))/3))))

    return '{:.1f}{}'.format(n / 10**(3 * millidx), millnames[millidx])

def getCol_Stats(runs, attrib):
    """Returns the column values over all of the samples"""
    ls = []
    for r in runs:
        ls.append((r['tumor']['id'], r['tumor'][attrib]))
        ls.append((r['normal']['id'], r['normal'][attrib]))
    vals = [x[1] for x in ls]
    mean = sum(vals)/len(vals)
    std = np.std(vals)
    return (ls, mean, std)

def printCell(sample, stats, attrib):
    n = sample[attrib]
    mean, std = stats[attrib]
    if n < (mean - 1.5*std): #1.5 ~ bottom 10%, 1 std ~ bottom 33%
        s = "red:%s" % millify(n)
    else:
        s = millify(n)
    return s

def main():
    usage = "USAGE: %prog -f [wes json file] -f [wes json file] ...  -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="mapping stats .csv file")
    optparser.add_option("-d", "--dir", help="section sub-dir")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.files or not options.dir or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #read in json data
    runs = [processJson(f) for f in options.files]

    #get the col values and calculate mean and std-dev for each attr
    stats = {}
    attrs = ['total_reads', 'mapped_reads', 'dedup_reads', 'mean_quality_score']
    for a in attrs:
        stats[a] = getCol_Stats(runs, a)

    #Write plot files:
    plots_dir = os.path.join(options.dir, 'plots')
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    #MAIN output.tsv
    out = open(options.output, "w")
    hdr = map(lambda x: prettyprint(x), attrs)
    hdr = ['&nbsp;'] + list(hdr)
    out.write("%s\n" % ",".join(hdr))
    row = []
    for (i, a) in enumerate(attrs):
        index = str(i).zfill(2)
        plot_name = "%s_%s_barh" % (index, a)
        csv_fname = "%s.csv" % plot_name
        det_fname = "%s_details.yaml" % (index)
        (ls, mean, std) = stats[a]
        csv = open(os.path.join(plots_dir, csv_fname), "w")
        hdr = ["Sample", prettyprint(a)]
        csv.write("%s\n" % ",".join(hdr))

        for (name, val) in ls:
            csv.write("%s\n" % ",".join([name,str(val)]))
        csv.close()
        #add link to image- clean path first --
        if options.dir.endswith("/"): #drop last /
            options.dir = options.dir[:-1]
        section = options.dir.split("/")[-1]
        plot_path = "img:%s/plots/%s.png" % (section, plot_name)
        if i == 0:
            out.write("&nbsp;,%s" % plot_path)
        else:
            out.write(",%s" % plot_path)
        #ALSO write details here!
    out.close()

if __name__ == '__main__':
    main()
