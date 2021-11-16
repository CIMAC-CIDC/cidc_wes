# Len Taing 2021 (TGBTG)
#MODULE: wes cohort report module
#RENAMES 'twist_maf_file' to 'filtered_maf_file'

import os
from yaml import dump as yaml_dump
from collections import OrderedDict

configfile: "config.cohort.yaml"

#cancerGeneList= "./cidc_wes/static/oncoKB/cancerGeneList.tsv"

def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for r in config['runs']:
        ls.append("jsons_filtered/%s.wes.json" % r)
    return ls

rule all:
    input:
        targets

def getJson(wildcards):
    json = config['runs'][wildcards.run][0] #should only be one
    return json

rule backport_sample_report:
    input:
        getJson
    params:
        sed_cmd="sed \'s/twist_maf_file/filtered_maf_file/g\'"
    output:
        "jsons_filtered/{run}.wes.json"
    shell:
        "cat {input} | {params.sed_cmd} > {output}"
