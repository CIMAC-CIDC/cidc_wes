#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import yaml

from string import Template

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

def addCondaPaths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    config['conda_root'] = conda_root
    config['wes_root'] = "%s/envs/wes" % conda_root

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
config = getRuns(config)
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
    #IMPORT all of the module targets
    ls.extend(align_targets(wildcards))
    ls.extend(metrics_targets(wildcards))
    ls.extend(recalibration_targets(wildcards))
    ls.extend(somatic_targets(wildcards))
    ls.extend(germline_targets(wildcards))
    ls.extend(checkmatch_targets(wildcards))
    ls.extend(coveragemetrics_targets(wildcards))
    ls.extend(copynumber_targets(wildcards))
    ls.extend(purity_targets(wildcards))
    ls.extend(neoantigen_targets(wildcards))
    ls.extend(optitype_targets(wildcards))
    #ls.extend(clonal_trial_targets(wildcards))
    #ls.extend(clonality_targets(wildcards))
    #ls.extend(report_targets(wildcards))
    return ls

rule target:
    input: 
        all_targets,

    message: "Compiling all output"
    

include: "./modules/align.snakefile"     # common align rules
include: "./modules/metrics.snakefile"   # ...
include: "./modules/recalibration.snakefile"      # ...
include: "./modules/somatic.snakefile" # ...
include: "./modules/germline.snakefile" # ...
include: "./modules/germlinematch.snakefile"
include: "./modules/coverage.snakefile" # ...
include: "./modules/copynumber.snakefile" # ...
include: "./modules/purity.snakefile" #...
include: "./modules/clonality.snakefile" # ...
include: "./modules/optitype.snakefile" #...
include: "./modules/neoantigen.snakefile"
#include: "./modules/report.snakefile"          # report module
