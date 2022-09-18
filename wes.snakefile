#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import json
import yaml

from string import Template

#LEN: ARE these center_targets now deprecated b/c we have twist targets?
#Dictionary of center targets
center_targets={'mocha':"./ref_files/hg38/target_beds/mocha.liftover.hg38.noContigs.bed",
                "mda": "./ref_files/hg38/target_beds/MDA.liftover.hg38.noContigs.bed",
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
    conda_root=os.environ.get("CONDA_ROOT")
    if not conda_root: #CONDA_ROOT not defined in envars, trying cmd line
        conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
    config['conda_root'] = conda_root
    #config['wes_root'] = "%s/envs/wes" % conda_root # DEPRECATED
    config['optitype_root'] = "%s/envs/optitype" % conda_root
    config['xhla_root'] = "%s/envs/xHLA" % conda_root
    config['sequenza_root'] = "%s/envs/sequenza" % conda_root
    config['pyclone_root'] = "%s/envs/pyclone" % conda_root
    config['vcf_root'] = "%s/envs/vcf" % conda_root


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

# we need skipped_modules to exist even if its empty
if 'skipped_modules' not in config:
    config['skipped_modules'] = []
    
# These modules can't be run in tumor only samples so they should be skipped in that case
if config.get('tumor_only'):
    for module in ['germline', 'purity', 'clonality', 'copynumber']:
        if module not in config['skipped_modules']:
            config['skipped_modules'].append(module)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)

#Finally check for 'remote_path' and 'expressions'
if 'remote_path' not in config:
    config['remote_path'] = ""
if 'expression_files' not in config:
    config['expression_files'] = {}
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
    if not config.get('tumor_only') and 'germline' not in config['skipped_modules']:
        ls.extend(germline_targets(wildcards))
    ls.extend(somatic_targets(wildcards))
    ls.extend(rna_targets(wildcards))
    return ls

def level2_sans_report(wildcards):
    skippable_module_dict = {
        'purity': purity_targets(wildcards),
        'clonality': clonality_targets(wildcards),
        'cnvkit': cnvkit_targets(wildcards),
        'neoantigen': neoantigen_targets(wildcards),
        'copynumber': copynumber_targets(wildcards),
        'msisensor2': msisensor2_targets(wildcards),
        'tcellextrect': tcellextrect_targets(wildcards),
    }
    # add optional modules to targets 
    ls = []
    for module in skippable_module_dict:
        if module not in config['skipped_modules']:
            ls.extend(skippable_module_dict[module])

            
    # add mandatory and special case modules to targets
    #ls.extend(copynumber_targets(wildcards))
    ls.extend(coveragemetrics_targets(wildcards))
    ls.extend(optitype_targets(wildcards))
    #if config.get('neoantigen_run_classII', False) and 'neoantigen' not in config['skipped_modules']:
    #Should run even if neoantigen is skipped
    if config.get('neoantigen_run_classII', False):
            ls.extend(xhla_targets(wildcards))
    return ls

def level2_targets(wildcards):
    ls = []
    ls.extend(level2_sans_report(wildcards))
    ls.extend(report_targets(wildcards))
    return ls

rule target:
    input:
        targets=all_targets,
        benchmarks="benchmarks.tar.gz"
    message: "Compiling all outputs"

rule level1:
    input: level1_targets
    message: "Compiling all LEVEL1 outputs"
    benchmark: "benchmarks/wes_level1_targets.txt"

rule level2:
    input: level2_targets
    message: "Compiling all LEVEL2 outputs"
    benchmark: "benchmarks/wes_level2_targets.txt"

rule tar_benchmarks:
    input: all_targets
    output: "benchmarks.tar.gz"
    shell: "tar -c benchmarks | gzip > {output}"


include: "./modules/align.snakefile"     # common align rules
include: "./modules/metrics.snakefile"   # sentieon basic qc metrics module
include: "./modules/recalibration.snakefile" # sentieon BQSR 
include: "./modules/somatic.snakefile" # somatic calling- General
include: "./modules/somatic_tnhaplotyper2.snakefile" # tnhaplotyper2 somatic caller
include: "./modules/somatic_tnsnv.snakefile" # tnsnv somatic caller
include: "./modules/somatic_tnscope.snakefile" # tnscope somatic caller (Default)
include: "./modules/germline.snakefile" # germline variant calling
include: "./modules/coverage.snakefile" # coverage metrics module
include: "./modules/copynumber.snakefile" # consensus CNV module
include: "./modules/purity.snakefile" # FACETS
include: "./modules/clonality.snakefile" # Sequenza and Pyclone6
include: "./modules/cnvkit.snakefile" # CNVkit
include: "./modules/optitype.snakefile" # HLA- Optitype
include: "./modules/xhla.snakefile" # HLA- xHLA
include: "./modules/hlahd.snakefile" # HLA- HLA-HD
include: "./modules/neoantigen.snakefile" #pvacseq neoantigen
include: "./modules/msisensor2.snakefile" #microsatellite instability
include: "./modules/tcellextrect.snakefile" #tcell estimation
include: "./modules/rna.snakefile" #somatic variants from RNA expression
include: "./modules/report.snakefile" # report module
#include: "./modules/report.cohort.snakefile" # cohort report module
