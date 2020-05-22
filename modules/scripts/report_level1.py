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

import math

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


def getAlignmentInfo(config):
    """Genereate the dictionary for the alignment section"""
    mapping_csv = "analysis/align/mapping.csv"
    ret = {'mapping_plot': 'wes_images/align/mapping.png',
           'mapping_csv': mapping_csv,
           'samples': []}
    for sample in config['samples']:
        sorted_bam = "analysis/align/%s/%s_sorted.bam" % (sample, sample)
        sorted_dedup_bam = "analysis/align/%s/%s_sorted.dedup.bam" % (sample, sample)
        gc_bias_file = " analysis/metrics/%s/%s_gc_metrics.txt" % (sample, sample)
        quality_score_file = " analysis/metrics/%s/%s_qd_metrics.txt" % (sample, sample)
        quality_by_cycle_file = " analysis/metrics/%s/%s_mq_metrics.txt" % (sample, sample)
        insert_size_file = " analysis/metrics/%s/%s_is_metrics.txt" % (sample, sample)
        tmp = {'name': sample,
               'gc_bias': 'wes_images/align/%s/%s_gcBias.png' % (sample, sample),
               
               'quality_score': 'wes_images/align/%s/%s_qualityScore.png' % (sample, sample), 
               'quality_by_cycle': 'wes_images/align/%s/%s_qualityByCycle.png' % (sample,sample),
               'insert_size': 'wes_images/align/%s/%s_insertSize.png' % (sample,sample),
               #FILES: are in the following form- (filename, filepath)
               'sorted_bam': (getFileName(sorted_bam),sorted_bam),
               'sorted_dedup_bam': (getFileName(sorted_dedup_bam),sorted_dedup_bam),
               'gc_bias_file': gc_bias_file,
               'quality_score_file': quality_score_file,
               'quality_by_cycle_file':quality_by_cycle_file,
               'insert_size_file': insert_size_file,
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

        #print out total as quantities of Millions of reads
        if int(d['total']) >= 1000:
            d['total'] = millify(int(d['total']))

        #GET the %_bases_above_50 from
        #analysis/metrics/{sample}/{sample}.{center}.mosdepth.region.summary.txt
        if 'cimac_center' in config and config['cimac_center']:
            center = config['cimac_center']
        else:
            #default broad
            center = 'broad'
        
        #coverage_file = "analysis/metrics/%s/%s_coverage_metrics.txt" % (sample,sample)
        coverage_file = "analysis/metrics/%s/%s.%s.mosdepth.region.summary.txt" % (sample,sample, center)
        #parse the coverage_file to get >50X
        if os.path.exists(coverage_file):
            cov_f = open(coverage_file)
            foo_hdr = cov_f.readline().strip().split(",")
            foo_vals = cov_f.readline().strip().split(",")
            cov_f.close()
            #NOTE: >50x is the second col
            bases_50x = "%.2f" % (float(foo_vals[1])*100.0)
            d['percent_bases_over_50'] = "%s %%" % bases_50x
        
        d['coverage_file'] = (getFileName(coverage_file),coverage_file)
        ret.append(d)
    #print(ret)
    f.close()
    return ret

def getSomaticSummaryTable(config):
    """Generates the mutation summary table"""
    #PARSE out the summary
    somatic_summary_table = {}
    mut_summary_f = "analysis/somatic/somatic_mutation_summaries.%s.csv" % config['somatic_caller']
    f = open(mut_summary_f)
    hdr = f.readline().strip().split(",")
    for l in f:
        tmp = l.strip().split(",")
        run_name = tmp[0]
        somatic_summary_table[run_name] = dict(zip(hdr,tmp))
    f.close()
    #print(somatic_summary_table)
    return somatic_summary_table

def getSomaticNonSynTable(config):
    """Generates the nonsynonymous summary table"""
    #PARSE out the summary
    somatic_nonsyn_table = {}
    mut_summary_f = "analysis/somatic/somatic_nonSynonymous_summaries.%s.csv" % config['somatic_caller']
    f = open(mut_summary_f)
    hdr = f.readline().strip().split(",")
    for l in f:
        tmp = l.strip().split(",")
        run_name = tmp[0]
        somatic_nonsyn_table[run_name] = dict(zip(hdr,tmp))
    f.close()
    return somatic_nonsyn_table

def getSNVSummaryTable(config):
    """Generates the Ref/Alt summary table, e.g. A>C, A>G counts"""
    #PARSE out the summary SNV
    somatic_SNV_table = {}
    for run in config['runs']:
        snv_summary_f = "analysis/somatic/%s/%s_%s_somatic_SNV_summaries.csv" % (run, run, config['somatic_caller'])
        somatic_SNV_table[run] = []
        f = open(snv_summary_f)
        hdr = f.readline().strip().split(",")
        #print(hdr)
        for l in f:
            tmp = l.strip().split(",")
            #print(tmp)
            somatic_SNV_table[run].append(tmp)
        f.close()
    #print(somatic_SNV_table)
    return somatic_SNV_table

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
               #NOTE: shirley just wants this called "{run}_{caller}.vcf}
               'filter_file': (getFileName(filter_file), "%s_%s.vcf" % (run, caller)),
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
                ('wes_level2.html','WES_Level2')]

    #SIDE-BAR
    alignment_sub = ['Mapping_Stats', 'GC_Bias_Plots',
                     'Quality_Score', 'Quality_by_Cycle']
    sidebar = [("alignment", "Alignment", alignment_sub),
               ('coverage','Coverage', []),
               ('germline',"Germline", []),
               ("somatic","Somatic", []),
    ]
    
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
    wes_report_vals['somatic_table'] = getSomaticSummaryTable(config)
    wes_report_vals['somatic_nonSyn_table'] = getSomaticNonSynTable(config)
    wes_report_vals['somatic_SNV_table'] = getSNVSummaryTable(config)
    wes_report_vals['somatic_summary_table_file'] = "analysis/somatic/somatic_mutation_summaries.%s.csv" % config['somatic_caller']

    #GERMLINE
    wes_report_vals['germline'] = getGermlineInfo(config)

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
