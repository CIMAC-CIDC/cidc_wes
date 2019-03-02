#module: clonality by Pyclone and vcftotsb by vcftools
#import os
#from string import Template


def clonalitypreprocess_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/align/%s/%s_recalibrated.bam" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def getNormal_sample(wildcards):
    return clonalitypreprocess_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return clonalitypreprocess_runsHelper(wildcards, 1)

def clonalitypreprocess_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s.seqz.txt.gz" % (run,run))
        ls.append("analysis/clonality/%s/%s.bin50.seqz.txt.gz" % (run,run))
        return ls

rule clonality_pre_all:
    input:
        clonalitypreprocess_targets

rule sequenza_bam2seqz:
    input:
        tumor_bam=getTumor_sample,
        normal_bam=getNormal_sample,
    output:
        sequenza_out="analysis/clonality/{run}/{run}.seqz.txt.gz"
    params:
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
        gc_file=config['gc_file'],
        ref=config['genome_fasta'],
        #chroms= ['chr%' % c for c in list(range(1,22+1))]
        #sequenz_path="/home/aashna/.local/bin", #LEN-add this 
    conda:
        "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_preclonality.txt"
    shell:
      "sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {output.sequenza_out}"

# rule sequenza_processing:
#     input:
#         tumor_bam=getTumor_sample,
#         normal_bam=getNormal_sample,
#     output:
#         sequenza_out="analysis/clonality/{run}/{run}.seqz.txt.gz"
#     params:
#         #JUST sample names - can also use the helper fns, e.g.
#         #normal = lambda wildcards: config['runs'][wildcards.run][0],
#         #tumor = lambda wildcards: config['runs'][wildcards.run][1],
#         gc_file=config['gc_file'],
#         ref=config['genome_fasta'],
#         #sequenz_path="/home/aashna/.local/bin", #LEN-add this 
#     conda:
#         "../envs/sequenza.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}_{run}_preclonality.txt"
#     shell:
#       "sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {output.sequenza_out}"

rule sequenza_preprocessing:
    input:
        completeseq="analysis/clonality/{run}/{run}.seqz.txt.gz",
    output:
        sequenzafinal_out="analysis/clonality/{run}/{run}.bin50.seqz.txt.gz"
    params:
        gc_file=config['gc_file'],
        ref=config['genome_fasta'],
        #sequenz_path="/home/aashna/.local/bin", #LEN-
    conda:
        "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_prestep2clonality.txt"
    shell:
      "sequenza-utils  seqz_binning --seqz {input.completeseq}  --window 50  -o {output.sequenzafinal_out}"



    
