#!/usr/bin/env python

"""Given a model from macs2 predictd, calculates the estimated fragment size
and the avg_sd

NOTE: these fns are directly translated from 
chilin2/modules/macs2_fragment/qc.py
"""
import os
import sys
from optparse import OptionParser

class SampleStats:
    """A class to store a set of statistics for the summary table
    NOTE: this information is aggregated over several .csv
    """
    def __init__(self, sampleID):
        self.sampleID = sampleID
    
def parseCSV(csv_file):
    """parses a CSV file, using the first line as a header (of column names)
    ASSUMES: sampleIDs are the first column
    returns a dictionary of dictionarys, e.g.
    {sampleID: {col1: val1, col2: val2, ...}}
    """
    ret = {}
    f = open(csv_file)
    cols = f.readline().strip().split(",")
    
    #read the rest
    for l in f:
        tmp = l.strip().split(',')
        #enter new record - don't enter the first col, which is sampleID
        ret[tmp[0]] = dict(zip(cols[1:], tmp[1:]))
    f.close()
    return ret

def addStat(d, csv, fieldToGet, fieldToStore):
    """ADD fieldToGet (value from csv) to d AS fieldToStore
    this fn has side-effects"""
    for s in list(d.keys()):
        if s in csv:
            d[s][fieldToStore] = csv[s][fieldToGet]
        else:
            d[s][fieldToStore] = "NA"
    

def main():
    usage = "USAGE: %prog -f [fastqc.csv] -m [mapping.csv] -p [pbc.csv] -r [fragSizes.csv] -b [bamRegionStats.csv]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--fastqc", help="fastqc.csv file")
    optparser.add_option("-m", "--mapping", help="mapping.csv file")
    optparser.add_option("-p", "--pbc", help="pbc.csv file")
    optparser.add_option("-r", "--fragSizes", help="fragSizes.csv file")
    optparser.add_option("-b", "--bamRegionStats", help="bamRegionStats.csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if (not options.fastqc or not options.mapping or not options.pbc \
            or not options.fragSizes or not options.bamRegionStats):
        optparser.print_help()
        sys.exit(-1)

    #print(",".join(["Sample","EstFragment","StdDev"]))

    #handle fastqc - get MedianQuality, store as FastQC
    tmp = parseCSV(options.fastqc)
    samples = list(tmp.keys())

    #initialize stats table--set keys to samples
    stats = {sampleID:{} for sampleID in samples}
    
    #get MedianQuality and set it as FastQC
    addStat(stats, tmp, 'MedianQuality', 'FastQC')

    #HANDLE mapping.csv
    tmp = parseCSV(options.mapping)
    addStat(stats, tmp, 'Total', 'TotalReads')
    addStat(stats, tmp, 'Mapped', 'MappedReads')
    addStat(stats, tmp, 'UniquelyMapped', 'UniqMappedReads')

    #HANDLE pbc.csv
    tmp = parseCSV(options.pbc)
    addStat(stats, tmp, 'Nd', 'UniqLoc4M')
    addStat(stats, tmp, 'N1', 'UniqLoc1read4M')
    #PBC = N1/Nd
    for s in samples:
        stats[s]['PBC'] = \
            int(stats[s]['UniqLoc1read4M']) / float(stats[s]['UniqLoc4M'])*100
        #FORMAT: convert to 2-decimal percentage
        stats[s]['PBC'] = "%.2f" % stats[s]['PBC']

    #HANDLE fragSizes.csv
    tmp = parseCSV(options.fragSizes)
    addStat(stats, tmp, 'MedianFrag', 'Fragment')

    #handle bamRegionStats.csv
    tmp = parseCSV(options.bamRegionStats)
    #add as percentage
    for s in samples:
        tot = int(tmp[s]['Total'])
        dhs = "%.2f" % (float(tmp[s]['DHS']) / tot *100)
        prom = "%.2f" % (float(tmp[s]['Promoter']) / tot *100)
        exon = "%.2f" % (float(tmp[s]['Exon']) / tot *100)
        stats[s]['DHS-Promoter-Exon4M'] = "%s/%s/%s" % (dhs, prom, exon)

    #print(stats)

    #OUTPUT- fields defines the column order
    fields = ['FastQC', 'TotalReads', 'MappedReads', 'UniqMappedReads',
              'UniqLoc4M', 'UniqLoc1read4M', 'PBC', 'DHS-Promoter-Exon4M']
    print(",".join(['Sample'] + fields))

    for s in samples:
        vals = [stats[s][f] for f in fields]
        print(",".join([s] + vals))

if __name__=='__main__':
    main()
