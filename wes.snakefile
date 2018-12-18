#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import yaml

from string import Template


def addPy2Paths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'chips_py2', 'lib', 'python2.7', 'site-packages')
    
    if not "python2" in config or not config["python2"]:
        config["python2"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'python2.7')

    if not "mdseqpos_path" in config or not config["mdseqpos_path"]:
        config["mdseqpos_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'MDSeqPos.py')

    if not "macs2_path" in config or not config["macs2_path"]:
        config["macs2_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'macs2')

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
#config = getRuns(config)
addPy2Paths_Config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)
#-----------------------------------------


#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------

def all_targets(wildcards):
    ls = []
    #IMPORT all of the module targets
    ls.extend(align_targets(wildcards))
    ls.extend(fastqc_targets(wildcards))
        
    #ls.extend(report_targets(wildcards))
    return ls

rule target:
    input: 
        all_targets,

    message: "Compiling all output"
    
include: "./modules/align_bwa.snakefile"     # rules specific to BWA
include: "./modules/align_common.snakefile"  # common align rules
include: "./modules/fastqc.snakefile"        # fastqc (sequence qual) rules
#include: "./modules/report.snakefile"        # report module
