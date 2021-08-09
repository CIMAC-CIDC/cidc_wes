# Len Taing 2021 (TGBTG)
#MODULE: wes cohort report module

import os
from yaml import dump as yaml_dump
from collections import OrderedDict

configfile: "config.cohort.yaml"

#cancerGeneList= "./cidc_wes/static/oncoKB/cancerGeneList.tsv"

def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for r in config['runs']:
        ls.append("jsons_addNeoantigen/%s.wes.json" % r)
    return ls

rule all:
    input:
        targets

def getJson(wildcards):
    json = config['runs'][wildcards.run][0] #should only be one
    neoantigen_results = config['neoantigen'][wildcards.run][0]
    tmp = {'json': json, 'neoantigen': neoantigen_results}
    #return json
    return tmp

rule backport_sample_report:
    input:
        unpack(getJson)
    #params:
        #KEY- retuires the center-specific targeted bed regions file
        #geneList = cancerGeneList,
    output:
        "jsons_addNeoantigen/{run}.wes.json"
    shell:
        "./cidc_wes/shell_scripts/report_backport/04_addNeoantigen.py -j {input.json} -f {input.neoantigen} -o {output}"
