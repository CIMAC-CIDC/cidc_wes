#MODULE: Indel and Base Realigner  by Sentieon
#import os
#from string import Template

_realigner_threads=15


def recal_runsHelper(wildcards, iindex, input_template):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    and a Python format string (ref: https://www.programiz.com/python-programming/methods/string/format)
    returns the template string with the run name"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append(input_template.format(sample=sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp
    
def recal_getNormal(wildcards):
    return recal_runsHelper(wildcards, 0, "analysis/align/{sample}/{sample}.realigned.bam")

def recal_getTumor(wildcards):
    return recal_runsHelper(wildcards, 1, "analysis/align/{sample}/{sample}.realigned.bam")

def recal_getNormal_table(wildcards):
    return recal_runsHelper(wildcards, 0, "analysis/align/{sample}/{sample}_prerecal_data.table")

def recal_getTumor_table(wildcards):
    return recal_runsHelper(wildcards, 1, "analysis/align/{sample}/{sample}_prerecal_data.table")


def recalibration_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/align/%s/%s.realigned.bam" % (sample,sample))
    	ls.append("analysis/align/%s/%s_prerecal_data.table" % (sample,sample))
    	ls.append("analysis/align/%s/%s_recalibrated.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_postrecal_data.table" % (sample,sample))
        ls.append("analysis/align/%s/%s_recal.csv" % (sample,sample))
        #ls.append("analysis/align/%s/%s_recal_plots.pdf" % (sample,sample))
        #ls.append("analysis/align/%s/%s.sort_recalibrated.bam" % (sample,sample))
    for run in config['runs']:
        ls.append("analysis/corealignments/%s/%s_tn_corealigned.bam" % (run,run))
    return ls

rule recalibration_all:
    input:
        recalibration_targets
    
rule Indel_realigner_sentieon:
    """indel realigner for uniquely  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.dedup.bam",
         bai="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
    output:
         realignbam="analysis/align/{sample}/{sample}.realigned.bam",
         realignbai="analysis/align/{sample}/{sample}.realigned.bam.bai"
    message:
         "INDEL REALIGNER: indel realigner for  mapped reads"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        #dbsnp=config['dbsnp'], #not used!
        mills=config['Mills_indels'],
        g1000=config['G1000_indels'],
    group: "recalibration"
    threads: _realigner_threads
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Indel_realigner_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} --algo Realigner -k {params.mills} -k {params.g1000} {output.realignbam}"""

rule Base_recalibration_precal_sentieon:
    """base recalibration for realigned files"""
    input:
        realignbam="analysis/align/{sample}/{sample}.realigned.bam",
    output:
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table",
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        recalibratedbai="analysis/align/{sample}/{sample}_recalibrated.bam.bai"
    message:
        " PRE BASE RECALIBRATION: base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: _realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_precal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.realignbam} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.prerecaltable} --algo ReadWriter {output.recalibratedbam}"""

rule Base_recalibration_postcal_sentieon:
    """post recalibration for realigned files"""
    input:
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        recalibrated_bai="analysis/align/{sample}/{sample}_recalibrated.bam.bai",
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table"
    output:
        postrecaltable="analysis/align/{sample}/{sample}_postrecal_data.table"
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: _realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_postcal_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.recalibratedbam} -q {input.prerecaltable} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.postrecaltable}"""

rule Base_recalibration_sentieon:
    """ recalibration for realigned files"""
    input:
        postrecaltable="analysis/align/{sample}/{sample}_postrecal_data.table",
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table"
    output:
        recalfile="analysis/align/{sample}/{sample}_recal.csv"
    message:
        "DIFF BASE RECALIBRATION: Difference in pre and post processing of realigned files"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
    threads: _realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{sample}/{sample}.Base_recalibration_sentieon.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} --algo QualCal --plot --before {input.prerecaltable} --after {input.postrecaltable} {output.recalfile}"""

rule corealignment:
    input:
        normal = recal_getNormal, 
        tumor = recal_getTumor,
        norm_recal = recal_getNormal_table,
        tumor_recal = recal_getTumor_table
    output:
        bam="analysis/corealignments/{run}/{run}_tn_corealigned.bam",
        bai="analysis/corealignments/{run}/{run}_tn_corealigned.bam.bai",
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        #dbsnp= config['dbsnp'], #not used
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: _realigner_threads
    group: "recalibration"
    benchmark:
        "benchmarks/recalibration/{run}/{run}.corealignment.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.tumor} -i {input.normal} -q {input.tumor_recal} -q {input.norm_recal} --algo Realigner -k {params.mills} -k {params.g1000} {output.bam}"""

