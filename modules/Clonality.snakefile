#module: clonality by Pyclone and vcftotsb by vcftools
#import os
#from string import Template


def clonality_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/somaticVariants/%s/%s_recalibrated.bam" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def getNormal_sample(wildcards):
    return clonality_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return clonality_runsHelper(wildcards, 1)

def clonality_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s_pyclone.tsv" % (run,run))
        return ls

rule clonality_all:
    input:
        clonality_targets

rule vcftotab_conversion:
    input:
        tumor_vcf=getTumor_sample
    output:
        tsvout="analysis/somaticVariants/{run}/{run}_{caller}.output.tsv"
    params:
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
         vcfperl="/mnt/ssd/wes/vcft/bin/vcf-to-tab"  #LEN :change  this 
    benchmark:
        "benchmarks/clonality/{run}/{run}_{caller}_clonality.txt"
    shell:
      """zcat {input.tumor_vcf} | {params.vcfperl} > {output.tsvout}"""


rule pyclone_clonality:
   input:
       intsv="analysis/somaticVariants/{run}/{run}_{caller}.output.tsv"
   output:
       pyclone_configyaml="analysis/clonality/{run}/{run}_{caller}_pyclone_analysis/config.yaml"
       pyclone_plots="analysis/clonality/{run}/{run}_{caller}_pyclone_analysis/plots/{run}.pdf"
       pyclone_tables="analysis/clonality/{run}/{run}_{caller}_pyclone_analysis/tables/{run}.tsv"
       pyclone_trace="analysis/clonality/{run}/{run}_{caller}_pyclone_analysis/trace/{run}.txt"
       pyclone_yaml="analysis/clonality/{run}/{run}_{caller}_pyclone_analysis/yaml/{run}.yaml"
  params:
       #JUST sample names - can also use the helper fns, e.g.
       #normal = lambda wildcards: getNormal_sample(wildcards)
       #normal = lambda wildcards: config['runs'][wildcards.run][0],
       tumor = lambda wildcards: config['runs'][wildcards.run][1],
       out_dir=lambda wildcards: "analysis/clonality/%s/%s_%s " % (wildcards.run, wildcards.run, wildcards.caller)

  benchmark:
       "benchmarks/clonality/{run}/{run}.pyclone.analysis.txt"
   shell:
       """PyClone run_analysis_pipeline --in_files {input.intsv} --working_dir {params.out_dir}"""
