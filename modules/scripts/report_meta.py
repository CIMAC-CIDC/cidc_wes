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

def getMetaInfo(config, wes_version_file):
    """Gets and populates a dictionary with the values required for the 
    meta pg"""

    #READ WES version from wes_version file
    f = open(wes_version_file)
    wes_version = f.read().strip()
    f.close()
    
    #HARD_code
    ref_version = "ver1.0 (build date: 20190911)"

    #FROM CONFIG
    assembly_version = config['assembly']
    sentieon_version = config['sentieon_path'].split("/")[-2]
    somatic_caller = "%s (sentieon)" % config['somatic_caller']


    #LEN: HARD code section
    vep_version="ensemble-vep (91.3)"
    facets_version="facets (0.5.14)"
    optitype_version="optitype (1.3.2)"
    neoantigen_caller = "pvactools (1.3.7)" #config['neoantigen_callers'] #LEN: hard-code
    epitope_lengths = config['neoantigen_epitope_lengths']
    vcftools_version="vcftools (0.1.16)"
    bcftools_version="bcftools (1.9)"

    
    #GET wes current tag
    cmd = "snakemake -v".split(" ")
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    snakemake_version = "snakemake (%s)" % output.decode("utf-8").strip()

    tmp = {'wes_version' : wes_version,
           'ref_version' : ref_version,
           'assembly_version': assembly_version,
           'sentieon_version': sentieon_version,
           'somatic_caller': somatic_caller,
           'vep_version': vep_version,
           'facets_version':facets_version,
           'optitype_version':optitype_version,
           'neoantigen_caller': neoantigen_caller,
           'epitope_lengths': epitope_lengths,
           'vcftools_version': vcftools_version,
           'bcftools_version': bcftools_version,
           'snakemake_version': snakemake_version,
    }
    #print(tmp)
    return tmp

def main():
    usage = "USAGE: %prog -c [config.yaml file] -o [output html file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--config", help="config.yaml file")
    optparser.add_option("-v", "--wes_ver_file", help="file that stores the wes version")
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
    template = templateEnv.get_template("wes_meta.html")

    #STANDARD nav bar list
    nav_list = [('wes_meta.html', 'WES_META'),
                ('wes_level1.html','WES_Level1'),
                ('wes_level2.html','WES_Level2')]

    pg_name = "WES_META"
    sidebar = [("meta", "META", [])]
    
    wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                       'page_name': pg_name}

    #META- get the dict and add it to wes_report_vals
    meta = getMetaInfo(config, options.wes_ver_file)
    wes_report_vals['meta'] = meta

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
