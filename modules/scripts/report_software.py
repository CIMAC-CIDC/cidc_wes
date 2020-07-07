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

def parseYaml(yaml_file):
    """Parses a yaml file and returns a dictionary"""
    ret = {}
    with open(yaml_file, 'r') as stream:
        try:
            ret = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return ret

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


def getSoftwareInfo(config, wes_versions_file):
    """Gets and populates a dictionary with the values required for the 
    meta pg"""

    assembly_version = config['assembly']
    sentieon_version = wes_versions_file.get('sentieon_version',"N/A")
    somatic_caller = "%s (sentieon)" % config['somatic_caller']

    #Tools versions
    vep_version = wes_versions_file.get('vep_version',"N/A")
    facets_version = wes_versions_file.get('facets_version',"N/A")
    optitype_version = wes_versions_file.get('optitype_version',"N/A")
    pvactools_version = wes_versions_file.get('pvactools_version',"N/A")
    vcftools_version = wes_versions_file.get('vcftools_version',"N/A")
    bcftools_version = wes_versions_file.get('bcftools_version',"N/A")
    snakemake_version = wes_versions_file.get('snakemake_version',"N/A")

    tmp = [("Software", "Version"),
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
    config=parseYaml(options.config)
    
    #PARSE wes version files
    wes_ver_file = parseYaml(options.wes_ver_file)

    #get wes versions as a list of tuples
    tmp = getSoftwareInfo(config, wes_ver_file)
    #write output
    out = open(options.output,'w')
    for (k,v) in tmp:
        out.write("%s\t%s\n" % (k,v))
    out.close()
    
if __name__ == '__main__':
    main()
