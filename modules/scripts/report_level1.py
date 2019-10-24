#!/usr/bin/env python3
"""Script to generate analysis/report/wes_level1.html"""

import os
import sys

import re
import subprocess
import yaml

from optparse import OptionParser
import jinja2
import pandas as pd

################### THIS fn is copied from wes.snakefile ######################
def getRuns(config):
    """parse metasheet for Run groupings"""
    ret = {}

    #LEN: Weird, but using pandas to handle the comments in the file
    #KEY: need skipinitialspace to make it fault tolerant to spaces!
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    #SKIP the hdr
    for l in f[1:]:
        tmp = l.strip().split(",")
        #print(tmp)
        ret[tmp[0]] = tmp[1:]

    #print(ret)
    config['runs'] = ret
    return config
###############################################################################
def getFileName(path):
    """Returns the file name of a given file path"""
    filename = path.split("/")[-1]
    return filename

def getAlignmentInfo(config):
    """Genereate the dictionary for the alignment section"""
    ret = {'mapping': 'wes_images/align/mapping.png',
           'samples': []}
    for sample in config['samples']:
        sorted_bam = "analysis/align/%s/%s_sorted.bam" % (sample, sample)
        sorted_dedup_bam = "analysis/align/%s/%s_sorted.dedup.bam" % (sample, sample)
        tmp = {'name': sample,
               #'mapping': 'wes_images/align/%s/mapping.png' % sample,
               'gc_bias': 'wes_images/align/%s/%s_gcBias.png' % (sample, sample),
               'quality_score': 'wes_images/align/%s/%s_qualityScore.png' % (sample, sample), 
               'quality_by_cycle': 'wes_images/align/%s/%s_qualityByCycle.png' % (sample,sample),
               'insert_size': 'wes_images/align/%s/%s_insertSize.png' % (sample,sample),
               #FILES: are in the following form- (filename, filepath)
               'sorted_bam': (getFileName(sorted_bam),sorted_bam),
               'sorted_dedup_bam': (getFileName(sorted_dedup_bam),sorted_dedup_bam),
        }
        ret['samples'].append(tmp)
    #print(ret)
    return ret

def getCoverageInfo(config, coverage_file):
    """Genereate the dictionary for the coverage section
    open file "analysis/metrics/all_sample_summaries.txt
    """
    ret = []
    f = open(coverage_file)
    hdr = f.readline().strip().split("\t")
    hdr[-1] = "percent_bases_over_50" #otherwise it's %_bases_above_50
    for l in f:
        tmp = l.strip().split("\t")
        d = dict(zip(hdr,tmp))
        #FILES--use sample_id as sample name
        sample = d['sample_id']
        coverage_file = "analysis/metrics/%s/%s_coverage_metrics.txt" % (sample,sample)
        d['coverage_file'] = (getFileName(coverage_file),coverage_file)
        ret.append(d)
    #print(ret)
    f.close()
    return ret

def getSomaticInfo(config):
    """Genereate the dictionary for the somatic section"""
    ret = []
    somatic_caller = config['somatic_caller']
    for run in config['runs']:
        #FILES:
        caller= config['somatic_caller']
        output_file="analysis/somatic/%s/%s_%s.output.vcf" % (run,run,caller)
        filter_file="analysis/somatic/%s/%s_%s.filter.vcf" % (run,run,caller)
        tmp = {'name': run,
               #NEED to simplify these names!!!!!!!
               'lego_plot': 'wes_images/somatic/%s/%s_%s.legoPlot.png' % (run, run, somatic_caller),
               'output_file': (getFileName(output_file), output_file),
               'filter_file': (getFileName(filter_file), filter_file),
        }
        ret.append(tmp)
    #print(tmp)
    return ret

def getGermlineInfo(config):
    """Genereate the dictionary for the germline section"""
    ret = []
    for run in config['runs']:
        #get the file to check for Tumor/Normal match
        #HARD-CODED relative path link which should be ok
        f = open("analysis/germline/%s/%s_vcfcompare.txt" % (run, run))
        lines=f.readlines()
        
        #STATISTIC is in 10th line
        line = lines[9].strip().split("\t") #GET 2nd to last and last elm
        #regex ref: https://stackoverflow.com/questions/4894069/regular-expression-to-return-text-between-parenthesis
        p1 = re.search(r'\((.*?)\)',line[-2]).group(1)
        p2 = re.search(r'\((.*?)\)',line[-1]).group(1)
        #NOTE: p1 and p2 are strils like XX.Y%--with the % trailing
        #So we need to remove the % character
        percent = (float(p1[:-1]) + float(p2[:-1]))/2.0
        #print(percent)
        tmp = {'name': run, 'percent': "%.2f" % percent}

        #FILES: get the _haplotyper.output.vcf for each of the samples
        hap_files=[]
        for sample in config['runs'][run]:
            haplotyper_file = "analysis/germline/%s/%s_haplotyper.output.vcf" % (sample, sample)
            hap_files.append([getFileName(haplotyper_file), haplotyper_file])
        #print(hap_files)
        tmp['haplotyper_files']= hap_files
        ret.append(tmp)
        f.close()
    #print(ret)
    return ret

def main():
    usage = "USAGE: %prog -c [config.yaml file] -o [output html file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--config", help="config.yaml file")
    optparser.add_option("-o", "--output", help="output html file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.config or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #PARSE config.yaml
    config={}
    with open(options.config, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    config = getRuns(config)
    
    #ASSUMING it's being run as WES project level
    templateLoader = jinja2.FileSystemLoader(searchpath="cidc_wes/report")
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template("wes_level1.html")

    pg_name = "WES_LEVEL_1"
    #STANDARD nav bar list
    nav_list = [('wes_meta.html', 'WES_META'),
                ('wes_level1.html','WES_Level1'),
                ('wes_level2.html','WES_Level2'),
                ('wes_level3.html','WES_Level3')]

    #SIDE-BAR
    alignment_sub = ['Mapping_Stats', 'GC_Bias_Plots',
                     'Quality_Score', 'Quality_by_Cycle']
    sidebar = [("alignment", "Alignment", alignment_sub),
               ('coverage','Coverage', []),
               ("somatic","Somatic", []),
               ('germline',"Germline", [])]
    
    wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                       'page_name': pg_name}

    #ALIGNMENT
    wes_report_vals['alignment'] = getAlignmentInfo(config)

    #COVERAGE
    #HARD-CODED relative path link which should work for now
    coverage_f = "analysis/metrics/all_sample_summaries.txt"
    wes_report_vals['coverage'] = getCoverageInfo(config, coverage_f)

    #SOMATIC
    wes_report_vals['somatic'] = getSomaticInfo(config)

    #GERMLINE
    wes_report_vals['germline'] = getGermlineInfo(config)

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
