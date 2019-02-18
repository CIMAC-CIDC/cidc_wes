#MODULE: Purity by Facets
#import os
_puritycalls_threads=8

def purity_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/align/%s/%s_recalibrated.bam" % (sample_name, sample_name))
        tmp.append("analysis/align/%s/%s_recalibrated.bam.bai" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    return tmp

def getNormal_sample(wildcards):
    return purity_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return purity_runsHelper(wildcards, 1)

def puritycalls_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/purity/%s/%s_purity_results.txt" % (sample,sample))
    return ls

rule puritycalls_all:
    input:
        puritycalls_targets

rule Puritycalls_Facets:
    """Get the  recalibrated bam files from  mapped reads"""
    input:
         NormalBam=getNormal_sample,
         TumorBam= getTumor_sample
    output:
         puritycalls="analysis/purity/{sample}/{sample}_purity_results.txt",
    message:
        "FACETS: purity calculations for bam file"
    params:
        index=config['genome_fasta'],
        index1=config['facets_vcftar']
    threads: _puritycalls_threads
    benchmark:
        "benchmarks/puritycalls/{sample}/{sample}.purityresults.txt"
    shell:
        """./snp-pileup -q15 -Q20  {params.index1} {output} {input.normalBam} {input.tumorBam}"""
