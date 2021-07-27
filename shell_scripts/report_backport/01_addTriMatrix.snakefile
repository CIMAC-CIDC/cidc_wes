# Len Taing 2021 (TGBTG)
#MODULE: wes cohort report module

import os
from yaml import dump as yaml_dump
from collections import OrderedDict

configfile: "config.cohort.yaml"
_bwa_index = './ref_files/hg38/bwa_gdc/GRCh38.d1.vd1.CIDC.fa'


def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for r in config['runs']:
        ls.append("json_backports/%s.wes.json" % r)
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
        #KEY- mutProfile needs the abs path to the reference
        ref = os.path.abspath(_bwa_index) 
    output:
        "json_backports/{run}.wes.json"
    shell:
        "./cidc_wes/shell_scripts/report_backport/01_addTriMatrix.py -j {input} -r {params.ref} -o {output}"
