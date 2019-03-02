#module: clonality by Pyclone and vcftotsb by vcftools
#import os
#from string import Template


def clonalitypost_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        #tmp.append("analysis/clonality/%s/%s.pyclone.tsv" % (sample_name, sample_name))
        tmp.append(sample_name)
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def getNormal_samplename(wildcards):
    return clonalitypost_runsHelper(wildcards, 0)

def getTumor_samplename(wildcards):
    return clonalitypost_runsHelper(wildcards, 1)

def clonality_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #ls.append("analysis/clonality/%s/%s_pyclone.tsv" % (run,run)),
        ls.append("analysis/clonality/%s/config.yaml" % (run,run))
        #ls.append("analysis/clonality/%s/plots/%s.pdf" % (run,run)),
        #ls.append("analysis/clonality/%s/tables/%s.tsv" % (run,run)),
        #ls.append("analysis/clonality/%s/trace/%s.txt" % (run,run)),
        #ls.append("analysis/clonality/%s/yaml/%s.yaml" % (run,run))
        #ls.append("analysis/clonality/pyclone/%s" % (run,run))
        return ls

rule clonality_all:
    input:
        clonality_targets


rule pyclone_finalprocessing:
    input:
        "analysis/clonality/{run}/{run}.pyclone.tsv"
    output:
        "analysis/clonality/{run}/config.yaml"
        #"analysis/clonality/{run}/plots/{run}.pdf",
        #"analysis/clonality/{run}/tables/{run}.tsv",
        #"analysis/clonality/{run}/trace/{run}.txt",
        #"analysis/clonality/{run}/yaml/{run}.yaml",
        #out_dir="analysis/clonality/pyclone"
    params:
        #out_dir=lambda wildcards: "analysis/clonality/%s/" % (wildcards.run, wildcards.run)
        out_dir="analysis/clonality"
        #pyclone_env="/home/taing/miniconda3/envs/pyclone/bin/"
    conda:
        "/mnt/ssd/wes/cidc_wes/envs/pyclone.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone.analysis.txt"
    shell:
        """PyClone run_analysis_pipeline --in_files {input}  --working_dir {params.out_dir}/{wildcards.run} """
