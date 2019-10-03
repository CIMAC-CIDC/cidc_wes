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

def getCNVInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    circos_plot = "wes_images/align/mapping.png" #LEN: CIRCOS plot here
    tmp = {"circos_plot": circos_plot}
    #print(tmp)
    return tmp

def getPurityInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    #read in the purity file for each RUN
    ret = []
    for run in config['runs']:
        f_name = "analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run)
        if os.path.isfile(f_name):
            f = open(f_name)
            hdr = f.readline().strip().split("\t")
            vals = dict(zip(hdr, f.readline().strip().split("\t")))
            #print(vals)
            tmp = {'name': run, 'purity': vals['purity']}
            ret.append(tmp)
        else:
            print("WARNING: expected %s, but it does not exist" % f_name)
    return ret

def getHLAInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    ret = []
    for sample in config['samples']:
        f_name = "analysis/optitype/%s/%s_result.tsv" % (sample,sample)
        if os.path.isfile(f_name):
            f = open(f_name)
            hdr = f.readline().strip().split("\t")
            vals = dict(zip(hdr, f.readline().strip().split("\t")))
            #print(vals)
            tmp = {'name': sample, 'HLA': vals}
            ret.append(tmp)
        else:
            print("WARNING: expected %s, but it does not exist" % f_name)
    return ret

def getNeoantigenInfo(config):
    """Gets and populates a dictionary with the values required for the page"""
    #read in the purity file for each RUN
    ret = []
    for run in config['runs']:
        neoantigen_plot = "wes_images/align/mapping.png" #LEN: plots here
        tmp = {'name': run, 'neoantigen_plot': neoantigen_plot}
        ret.append(tmp)
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
    template = templateEnv.get_template("wes_level2.html")

    #STANDARD nav bar list
    nav_list = [('wes_meta.html', 'WES_META'),
                ('wes_level1.html','WES_Level1'),
                ('wes_level2.html','WES_Level2'),
                ('wes_level3.html','WES_Level3')]

    #MODIFY this!
    pg_name = "WES_LEVEL_2"
    sidebar = [("cnv", "CNV", []),
               ("purity", "Purity", []),
               ('hla','HLA', []),
               ("neoantigen","Neoantigen", [])]
    
    wes_report_vals = {'top_nav_list':nav_list, 'sidebar_nav': sidebar,
                       'page_name': pg_name}

    #CNV- get the dict and add it to wes_report_vals
    cnv = getCNVInfo(config)
    wes_report_vals['cnv'] = cnv

    #PURITY
    purity = getPurityInfo(config)
    wes_report_vals['purity'] = purity

    #HLA
    hla = getHLAInfo(config)
    wes_report_vals['hla'] = hla

    #Neoantigen
    neoantigen = getNeoantigenInfo(config)
    wes_report_vals['neoantigen'] = neoantigen

    template.stream(wes_report_vals).dump(options.output)  
        
if __name__ == '__main__':
    main()
