#MODULE: Indel and Base Realigner  by Sentieon
#import os
#from string import Template

_realigner_threads=8
_dbsnp="/cluster/asahu/mutation_calling/MDAnderson/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
_Mills_indels="/cluster/asahu/mutation_calling/MDAnderson/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
_G1000_indels="/cluster/asahu/mutation_calling/script/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz"


def runsHelper(wildcards, iindex, input_template):
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
    

def getNormal(wildcards):
    return runsHelper(wildcards, 0, "analysis/align/{sample}/{sample}.realigned.bam")

def getTumor(wildcards):
    return runsHelper(wildcards, 1, "analysis/align/{sample}/{sample}.realigned.bam")

def getNormal_recal(wildcards):
    return runsHelper(wildcards, 0, "analysis/align/{sample}/{sample}_prerecal_data.table")

def getTumor_recal(wildcards):
    return runsHelper(wildcards, 1, "analysis/align/{sample}/{sample}_prerecal_data.table")


def recalibration_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #TODO- fill this in
    return ls

rule recalibration_all:
    input:
        recalibration_targets

rule Indel_realigner_sentieon:
    """indel realigner for uniquely  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.unique.dedup.sorted.bam"
    output:
         realignbam="analysis/align/{sample}/{sample}.realigned.bam"
    message:
         "INDEL REALIGNER: indel realigner for  mapped reads"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads: _realigner_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.bam} --algo Realigner -k {params.mills} -k {params.g1000} {output.realignbam}"""

rule Base_recalibration_precal_sentieon:
    """base recalibration for realigned files"""
    input:
        realignbam="analysis/align/{sample}/{sample}.realigned.bam",
    output:
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table",
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam"
    message:
        " PRE BASE RECALIBRATION: base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
	dbsnp=_dbsnp,
	mills=_Mills_indels,
	g1000=_G1000_indels,
    threads: _realigner_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.realignbam} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.prerecaltable} --algo Readwriter {output.recalibratedbam}"""

rule Base_recalibration_postcal_sentieon:
    """post recalibration for realigned files"""
    input:
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table"
    output:
        postrecaltable="analysis/align/{sample}/{sample}_recal_data.table"
    message:
        "POST BASE RECALIBRATION: post base recalibration for  realigned files"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
	dbsnp=_dbsnp,
	mills=_Mills_indels,
	g1000=_G1000_indels,
    threads: _realigner_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.recalibratedbam} -q {input.prerecaltable} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000}  {output.postrecaltable}"""

rule Base_recalibration_sentieon:
    """ recalibration for realigned files"""
    input:
        postrecaltable="analysis/align/{sample}/{sample}_recal_data.table",
        prerecaltable="analysis/align/{sample}/{sample}_prerecal_data.table"
    output:
        recalfile="analysis/align/{sample}/{sample}_recal.csv"
    message:
        "DIFF BASE RECALIBRATION: Difference in pre and post processing of realigned files"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
    threads: _realigner_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} --algo QualCal --plot --before {input.prerecaltable} --after {input.postrecaltable} {output.recalfile}"""

rule Base_recalibration_plot:
    """base realigner plotter for recalibrated files"""
    input:
        recalfile="analysis/align/{sample}/{sample}_recal.csv"
    output:
        recalplot="analysis/align/{sample}/{sample}_recal_plots.pdf"
    message:
        "BASE RECALIBRATED PLOTS: base realigner plotting"
    params:
        #index=config['genome_fasta'],
        index1=config['sentieon_path'],
    threads: _realigner_threads
    shell:
        """{params.index1}/sentieon plot bqsr -o {output.recalplot}  {input.recalfile}"""

rule corealignment:
    input:
        norm = getNormal, 
        tumor = getTumor,
        norm_recal = getNormal_recal,
        tumor_recal = getTumor_recal
    output:
        "analysis/runs/{run}/{run}_tn_corealigned.bam"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads: _realigner_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.tumor} -i {input.normal} -q {input.tumor_recal} -q {input.norm_recal} --algo Realigner -k {params.mills} -k {params.g1000} {output}"""
