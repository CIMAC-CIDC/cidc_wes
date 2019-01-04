#!/usr/bin/env python

"""Given a model from macs2 predictd, calculates the estimated fragment size
and the avg_sd

NOTE: these fns are directly translated from 
chilin2/modules/macs2_fragment/qc.py
"""
import os
import sys
from optparse import OptionParser

import math

def get_size(rscript):
    """PARSES the model R file that is generated by macs2, and instead 
    of making an R script, it translates it to a python dictionary: values"""
    values = {}
    with open(rscript) as model:
        for line in model:
            if line.startswith("p <- "):
                #CONVERT to python tuple
                values["positive"] = eval(line.replace("p <- c", ""))
            if line.startswith("m <- "):
                values["minus"] = eval(line.replace("m <- c", ""))
            if line.startswith("x <- "):
                #NOTE: this can be done outside the loop!
                p = values['positive']
                r = int((len(p) - 1)/2)
                values["x"] = range(-1*r, r)
            if line.startswith("xcorr"):
                values['xcorr'] = eval(line.replace("xcorr <- c", ""))
            if line.startswith("ycorr"):
                values['ycorr'] = eval(line.replace("ycorr <- c", ""))
    return values

def calcFrag(values):
    """Given a set of parsed values that were generated by macs2 predictd,
    see get_size fn, this fn calculates the estimated fragment size and 
    the sd.
    
    **IMPORTANT: this fn is a python translation of the R code in 
    chilin2/modules/macs2_fragment/qc.py -- stat_frag_std

    RETURNS: (estimated frag, sd)
    """
    #calculate the estimated frag size: xmax
    ymax = max(values['ycorr'])
    i_ymax = values['ycorr'].index(ymax)
    xmax = values['xcorr'][i_ymax]
    #print(ymax, xmax)

    #find expected
    p_expect=sum([x* p/100.0 for (x,p) in zip(values['x'],values['positive'])])
    m_expect=sum([x* m/100.0 for (x,m) in zip(values['x'],values['minus'])])
    #print(p_expect, m_expect)

    #calc sd
    p_sd = math.sqrt(sum([((x - p_expect)**2)* p/100.0 \
                              for (x,p) in zip(values['x'],values['positive'])]))
    m_sd = math.sqrt(sum([((x - m_expect)**2)* m/100.0 \
                              for (x,m) in zip(values['x'],values['minus'])]))
    #print(p_sd, m_sd)

    #FINAL avg std error
    avg_sd = (p_sd + m_sd) /2.0

    return (xmax, avg_sd)

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of fastqc/{sample}.csv files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    print(",".join(["Sample","MedianFrag","StdDev"]))

    for f in options.files:
        #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
        sampleID = f.strip().split("/")[-1].split('.')[0]
        if sampleID.endswith('_fragModel'):
            sampleID = sampleID.replace('_fragModel','')

        values = get_size(f)
        (frag, sd) = calcFrag(values)
        frag = "%.2f" % frag
        sd = "%.2f" % sd
        print(",".join([sampleID,frag,sd]))

if __name__=='__main__':
    main()