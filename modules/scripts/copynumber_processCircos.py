#!/usr/bin/env python
"""
smooths sentieon _cnvcalls.txt.tn.tsv values by certain window size

NOTE: only works for human samples!
BELOW are hg38 chrom sizes-
chr1	248956422
chr2	242193529
chr3	198295559
chr4	190214555
chr5	181538259
chr6	170805979
chr7	159345973
chr8	145138636
chr9	138394717
chr10	133797422
chr11	135086622
chr12	133275309
chr13	114364328
chr14	107043718
chr15	101991189
chr16	90338345
chr17	83257441
chr18	80373285
chr19	58617616
chr20	64444167
chr21	46709983
chr22	50818468
chrX	156040895
chrY	57227415
"""

import os
import sys
from optparse import OptionParser

            
def main():
    usage = "USAGE: %prog -f [tab delim txt] (-w window size)"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--infile", help="file to filter")
    optparser.add_option("-w", "--window", help="window size (500000)", default=500000)
    (options, args) = optparser.parse_args(sys.argv)

    if not options.infile:
        optparser.print_help()
        sys.exit(-1)

    window = int(options.window)
    chromLens = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
                 'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
                 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
                 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
                 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                 'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
                 'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
                 'chr22': 50818468 , 'chrX': 156040895, 'chrY': 57227415}
    order = ['chr1', 'chr2', 'chr3',  'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
             'chr9','chr10','chr11', 'chr12','chr13','chr14','chr15','chr16',
             'chr17','chr18', 'chr19',  'chr20','chr21','chr22','chrX','chrY']
    
    #STORE them as a dictionary (chroms) lists of lists
    #initialize the dictionary to empty arrays
    cnv = {}
    for chrom in chromLens:
        n = int(chromLens[chrom]/window) + 1
        cnv[chrom] = [[] for i in range(n)]

    #for chrom in cnv:
    #    print("%s:%s" % (chrom, len(cnv[chrom])))

    f = open(options.infile)

    hdr = f.readline().strip().split("\t")
    
    #rest of the file
    for l in f:
        tmp = l.strip().split("\t")
        chrom = tmp[0]
        start = tmp[1]
        end = tmp[2]
        val = float(tmp[4])
        
        #calculate which bin it should belong
        start_bin = int(int(start) / window)
        #print(start_bin)
        #end_bin = int(end /window)

        #NO double counting!--assign the cnv to the start_bin
        cnv[chrom][start_bin].append(val)

    #convert this into a cnv output
    #hdr
    print("\t".join([hdr[0],hdr[1],hdr[2],hdr[4]]))
    for chrom in order:
        for i in range(int(chromLens[chrom]/window) + 1):
            #convert to hs
            tmp = "hs%s" % chrom[3:]
            start = i*window
            end = start + window if (start + window < chromLens[chrom]) else chromLens[chrom]
            
            mean = sum(cnv[chrom][i])/len(cnv[chrom][i]) if cnv[chrom][i] else 0
            print("\t".join([tmp, str(start), str(end), "%.4f" % mean]))
    f.close()


if __name__=='__main__':
    main()
