# Len Taing 2021 (TGBTG)
#MODULE: wes cohort report module

import os
from yaml import dump as yaml_dump
from collections import OrderedDict

configfile: "config.cohort.yaml"

#Dictionary of center targets
center_targets={'mocha':"./ref_files/hg38/target_beds/mocha.liftover.hg38.noContigs.bed",
                "mda": "./ref_files/hg38/target_beds/MDA.liftover.hg38.noContigs.bed",
                "broad":"./ref_files/hg38/target_beds/broad.liftover.hg38.bed"}

def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for r in config['runs']:
        ls.append("jsons_addSomaticSummary/%s.wes.json" % r)
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
        targets = center_targets[config.get('center', 'broad')],
    output:
        "jsons_addSomaticSummary/{run}.wes.json"
    shell:
        "./cidc_wes/shell_scripts/report_backport/02_addSomaticSummary.py -j {input} -t {params.targets} -o {output}"
