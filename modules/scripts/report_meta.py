#!/usr/bin/env python3
"""Script to generate wes run information"""

import os
import sys

import re
import subprocess
import yaml

from optparse import OptionParser
#import jinja2
#import pandas as pd

def getWESCommit():
    """Tries to get the wes commit string by system calls"""
    
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
    #wes_version = "WES %s (commit: %s)" % (wes_tag, wes_commit)
    #CHANGE back to CWD
    os.chdir(wd)
    return (wes_tag, wes_commit)


def getMetaInfo(config, wes_versions_file):
    """Gets and populates a dictionary with the values required for the 
    meta pg"""

    #WES version
    (wes_tag, wes_commit) = getWESCommit()
    wes_version = "WES %s (commit: %s)" % (wes_tag, wes_commit)
    wes_ref_version = config.get('wes_ref_version', None)
    #if not defined, then get the default from wes_versions_file
    if not wes_ref_version:
        wes_ref_version = wes_versions_file.get('wes_ref_version','N/A')
        
    wes_image = config.get('wes_image', None)
    if not wes_image:
        wes_image = wes_versions_file.get('wes_image','N/A')

    #FROM CONFIG
    assembly_version = config['assembly']
    sentieon_version = config['sentieon_path'].split("/")[-2]
    somatic_caller = "%s (sentieon)" % config['somatic_caller']


    #Tools versions
    vep_version = wes_versions_file.get('vep_version',"N/A")
    facets_version = wes_versions_file.get('facets_version',"N/A")
    optitype_version = wes_versions_file.get('optitype_version',"N/A")
    pvactools_version = wes_versions_file.get('pvactools_version',"N/A")
    vcftools_version = wes_versions_file.get('vcftools_version',"N/A")
    bcftools_version = wes_versions_file.get('bcftools_version',"N/A")
    snakemake_version = wes_versions_file.get('snakemake_version',"N/A")

    tmp = [('WES Version', wes_version),
           ('WES Reference Files', wes_ref_version),
           ('WES Tools Image', wes_image),
           ('Assembly Version', assembly_version),
           ('Sentieon Version', sentieon_version),
           ('Somatic Caller', somatic_caller),
           ('Ensembl VEP Version', vep_version),
           ('Facets Version', facets_version),
           ('Optitype Version (HLA caller)', optitype_version),
           ('Pvactools Version (neoantigen caller)', pvactools_version),
           ('Vcftools Version', vcftools_version),
           ('Bcftools Version', bcftools_version),
           ('Snakemake Version', snakemake_version)]

    #print(tmp)
    return tmp

def main():
    usage = "USAGE: %prog -c [config.yaml file] -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--config", help="config.yaml file")
    optparser.add_option("-v", "--wes_ver_file", help="yaml file that stores wes software versions ")
    optparser.add_option("-o", "--output", help="output tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if not options.config or not options.wes_ver_file or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #PARSE config.yaml
    config={}
    with open(options.config, 'r') as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    #PARSE wes version files
    wes_ver_file = {}
    with open(options.wes_ver_file, 'r') as stream:
        try:
            wes_ver_file = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    #get wes versions as a list of tuples
    tmp = getMetaInfo(config, wes_ver_file)
    #write output
    out = open(options.output,'w')
    for (k,v) in tmp:
        out.write("%s\t%s\n" % (k,v))
    out.close()
    
if __name__ == '__main__':
    main()
