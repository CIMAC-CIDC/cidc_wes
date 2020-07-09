#!/usr/bin/env python3
"""Script to generate wes run information"""

import os
import sys

from optparse import OptionParser

def composeRow(sample, image_path, isTumor=True):
    """Given a run, sample name, and whether it's tumor or normal, generate 
    the row for the tsv"""
    #relative path is the path relative to report.html, i.e. removing
    #"analysis/report"
    rel_path = "/".join(image_path.split("/")[2:])
    ttype = "T" if isTumor else "N"
    tmp = "\t".join(["%s(%s)" % (sample, ttype),
                     "img:%s" % os.path.join(rel_path, "%s_gcBias.png" % sample),
                     "img:%s" % os.path.join(rel_path, "%s_qualityScore.png" % sample),
                     "img:%s" % os.path.join(rel_path, "%s_qualityByCycle.png" % sample),
                     "img:%s" % os.path.join(rel_path, "%s_insertSize.png" % sample)])
    #print(tmp)
    return tmp

def main():
    usage = "USAGE: %prog -n [name of normal sample] -t [name of tumor sample]  -i {path to where the image files are stored} -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-n", "--normal", help="name of normal sample")
    optparser.add_option("-t", "--tumor", help="name of tumor sample")
    optparser.add_option("-p", "--path", help="path to where the images are stored")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.normal or not options.tumor or not options.path or not options.output:
        optparser.print_help()
        sys.exit(-1)
    hdr = "\t".join(["Sample","GC Plot", "Quality Score", "Qualty By Cycle",
                     "Insert Size"])
    out = open(options.output, "w")
    out.write("%s\n" % hdr)
    out.write("%s\n" % composeRow(options.tumor, options.path, True))
    out.write("%s\n" % composeRow(options.normal, options.path, False))
    out.close()

if __name__ == '__main__':
    main()
