# Len Taing 2021 (TGBTG)
#MODULE: wes cohort report module

import os
from yaml import dump as yaml_dump
from collections import OrderedDict

configfile: "config.cohort.yaml"

cancerGeneList= "./cidc_wes/static/oncoKB/cancerGeneList.tsv"

def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for r in config['runs']:
        ls.append("jsons_addGeneList/%s.wes.json" % r)
    return ls

rule all:
    input:
        targets

def getJson(wildcards):
    s = config['runs'][wildcards.run]
    return s[0] #should only be one

rule backport_sample_report:
    input:
        getJson
    params:
        #KEY- retuires the center-specific targeted bed regions file
        geneList = cancerGeneList,
    output:
        "jsons_addGeneList/{run}.wes.json"
    shell:
        "./cidc_wes/shell_scripts/report_backport/03_addGeneList.py -j {input} -l {params.geneList} -o {output}"
