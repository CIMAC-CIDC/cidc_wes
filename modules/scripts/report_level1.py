#!/usr/bin/env python3
"""Script to generate analysis/report/wes_level1.html"""

import os
import sys
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

def getMetaInfo(config):
    """Gets and populates a dictionary with the values required for the 
    meta pg"""
    
    # tmp = {'wes_version' : "v1.1 (commit: d8c124c)",
    #        'ref_version' : "ver1.0 (build date: 20190911)",
    #        "assembly_version": "GDC hg38",
    #        "sentieon_version": "201808.05",
    #        "somatic_caller": "tnscope (sentieon)",
    #        "neoantigen_callers": "MHCflurry NetMHCcons MHCnuggetsII",
    #        "epitope_lengths": "8,9,10,11",
    #        "snakemake_version": "5.4.5"}

    #CHANGE to cidc_wes, but first store this
    wd = os.getcwd()
    os.chdir("cidc_wes")
    #GET wes commit string--first six letters of this cmd
    cmd = "git show --oneline -s".split(" ")
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    wes_commit = output[:6].decode("utf-8") 

    #GET wes current tag
    cmd = "git describe --tags".split(" ")
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    wes_tag = output.decode("utf-8").split("-")[0]
    wes_version = "%s (commit: %s)" % (wes_tag, wes_commit)
    #CHANGE back to CWD
    os.chdir(wd)

    #HARD_code
    ref_version = "ver1.0 (build date: 20190911)"

    #FROM CONFIG
    assembly_version = config['assembly']
    sentieon_version = config['sentieon_path'].split("/")[-2]
    somatic_caller = "%s (sentieon)" % config['somatic_caller']
    neoantigen_callers = config['neoantigen_callers']
    epitope_lengths = config['neoantigen_epitope_lengths']

    #GET wes current tag
    cmd = "snakemake -v".split(" ")
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    snakemake_version = output.decode("utf-8").strip()

    tmp = {'wes_version' : wes_version,
           'ref_version' : ref_version,
           'assembly_version': assembly_version,
           'sentieon_version': sentieon_version,
           'somatic_caller': somatic_caller,
           'neoantigen_callers': neoantigen_callers,
           'epitope_lengths': epitope_lengths,
           'snakemake_version': snakemake_version}
    #print(tmp)
    return tmp

def getAlignmentInfo(config):
    """Genereate the dictionary for the alignment section"""
    ret = []
    for sample in config['samples']:
        tmp = {'name': sample,
               'mapping': 'analysis/report/wes_images/align/%s/mapping.png' % sample,
               'gc_bias': 'analysis/report/wes_images/align/%s/%s_metrics_1.png' % (sample, sample),
               'quality_score': 'analysis/report/wes_images/align/%s/%s_metrics_2.png' % (sample, sample), 
               'quality_by_cycle': 'analysis/report/wes_images/align/%s/%s_metrics_3.png' % (sample,sample)
        }
        ret.append(tmp)
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
        ret.append(dict(zip(hdr,tmp)))
    #print(ret)
    return ret

def getSomaticInfo(config):
    """Genereate the dictionary for the somatic section"""
    ret = []
    somatic_caller = config['somatic_caller']
    for run in config['runs']:
        tmp = {'name': run,
               #NEED to simplify these names!!!!!!!
               'lego_img': 'analysis/report/wes_images/somatic/%s/%s_%s.output_1.png' % (run, run, somatic_caller)}
        ret.append(tmp)
    #print(tmp)
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
    nav_list = [('wes_level1.html','WES_Level1'),
                ('wes_level2.html','WES_Level2'),
                ('wes_level3.html','WES_Level3')]

    #SIDE-BAR
    alignment_sub = ['Mapping_Stats', 'GC_Bias_Plots',
                     'Quality_Score', 'Quality_by_Cycle']
    sidebar = [("meta", "META", []),
               ("alignment", "Alignment", alignment_sub),
               ('coverage','Coverage', []),
               ("somatic","Somatic", []),
               ('germline',"Germline", [])]
    
    wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                       'page_name': pg_name}

    #META- get the dict and add it to wes_report_vals
    meta = getMetaInfo(config)
    wes_report_vals['meta'] = meta

    #ALIGNMENT
    wes_report_vals['alignment'] = getAlignmentInfo(config)

    #COVERAGE
    #HARD-CODED relative path link which should work for now
    coverage_f = "analysis/metrics/all_sample_summaries.txt"
    wes_report_vals['coverage'] = getCoverageInfo(config, coverage_f)

    #SOMATIC
    wes_report_vals['somatic'] = getSomaticInfo(config)

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
