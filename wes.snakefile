#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import json
import yaml

from string import Template

#Dictionary of center targets
center_targets={'mocha':"./ref_files/hg38/target_beds/mocha.liftover.hg38.bed",
                "mda": "./ref_files/hg38/target_beds/MDA.liftover.hg38.bed",
                "broad":"./ref_files/hg38/target_beds/broad.liftover.hg38.bed"}

def getRunsCohorts(config):
    """parse metasheet for Run groupings"""
    ret = {}

    #LEN: Weird, but using pandas to handle the comments in the file
    #KEY: need skipinitialspace to make it fault tolerant to spaces!
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    #SKIP the hdr
    cohorts = {}
    for l in f[1:]:
        tmp = l.strip().split(",")
        #print(tmp)
        ret[tmp[0]] = tmp[1:3]

        if len(tmp) > 3:
            chort=tmp[3] #it's the 4th col
            if chort in cohorts:
                cohorts[chort].append(tmp[0]) #add run to the cohort set
            else:
                #new cohort
                cohorts[chort] = [tmp[0]]

    #Find any cohorts that are singletons and print warning
    delete_these = []
    for c in cohorts:
        if len(cohorts[c]) < 2:
            print("WES WARNING: cohort group %s ignored because it has too few runs in the set (at least 2 are required).  Please correct your metasheet" % c)
            delete_these.append(c)
    #REMOVE them
    for c in delete_these:
        del cohorts[c] #not valid cohort, remove

    #print(ret)
    config['runs'] = ret

    config['cohorts'] = cohorts if cohorts else None
    #print(config['cohorts'])
    return config

def addCondaPaths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    config['conda_root'] = conda_root
    config['wes_root'] = "%s/envs/wes" % conda_root
    config['optitype_root'] = "%s/envs/optitype" % conda_root
    config['xhla_root'] = "%s/envs/xHLA" % conda_root
    config['sequenza_root'] = "%s/envs/sequenza" % conda_root
    config['pyclone_root'] = "%s/envs/pyclone" % conda_root


def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
config = getRunsCohorts(config)
addCondaPaths_Config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)

#Finally check for 'remote_path'
if 'remote_path' not in config:
    config['remote_path'] = ""
#-----------------------------------------


#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------

def all_targets(wildcards):
    ls = []
    #using level helper fns--so I Don't repeat myself!
    lvl1_targets = level1_targets(wildcards)
    lvl2_targets = level2_targets(wildcards)
    ls.extend(lvl1_targets)
    ls.extend(lvl2_targets)
    #print("\n".join(ls))
    return ls

def level1_targets(wildcards):
    ls = []
    ls.extend(align_targets(wildcards))
    ls.extend(metrics_targets(wildcards))
    ls.extend(recalibration_targets(wildcards))
    ls.extend(germline_targets(wildcards))
    ls.extend(somatic_targets(wildcards))
    return ls

def level2_targets(wildcards):
    ls = []
    ls.extend(coveragemetrics_targets(wildcards))
    ls.extend(copynumber_targets(wildcards))
    ls.extend(purity_targets(wildcards))
    ls.extend(neoantigen_targets(wildcards))
    ls.extend(optitype_targets(wildcards))
    if 'neoantigen_run_classII' in config and config['neoantigen_run_classII']:
        ls.extend(xhla_targets(wildcards))
    ls.extend(clonality_targets(wildcards))
    ls.extend(report_targets(wildcards))
    return ls

rule target:
    input: 
        all_targets,
    message: "Compiling all outputs"
    benchmark: "benchmarks/all_wes_targets.txt"

rule level1:
    input: level1_targets
    message: "Compiling all LEVEL1 outputs"
    benchmark: "benchmarks/wes_level1_targets.txt"

rule level2:
    input: level2_targets
    message: "Compiling all LEVEL2 outputs"
    benchmark: "benchmarks/wes_level2_targets.txt"

include: "./modules/align.snakefile"     # common align rules
include: "./modules/metrics.snakefile"   # ...
include: "./modules/recalibration.snakefile"      # ...
include: "./modules/somatic.snakefile" # ...
include: "./modules/somatic_tnhaplotyper2.snakefile" # ...
include: "./modules/somatic_tnsnv.snakefile" # ...
include: "./modules/somatic_tnscope.snakefile" # ...
include: "./modules/germline.snakefile" # ...
include: "./modules/coverage.snakefile" # ...
include: "./modules/copynumber.snakefile" # ...
include: "./modules/purity.snakefile" #...
include: "./modules/clonality.snakefile" # ...
include: "./modules/optitype.snakefile" #...
include: "./modules/xhla.snakefile" #....
include: "./modules/neoantigen.snakefile"
include: "./modules/report.snakefile" # report module
