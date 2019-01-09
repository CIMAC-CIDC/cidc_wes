#MODULE: Somatic Variant calls by Sentieon
#import os
#from string import Template

_somaticcall_threads=8
_dbsnp="/cluster/asahu/mutation_calling/MDAnderson/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
_Mills_indels="/cluster/asahu/mutation_calling/MDAnderson/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
_G1000_indels="/cluster/asahu/mutation_calling/script/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz"


def runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    and a Python format string (ref: https://www.programiz.com/python-programming/methods/string/format)
    returns the template string with the run name"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
	#config['samples'][sample_name] = the fastq file paths defined in config.yaml
	tmp.extend(config['samples'][sample_name]) #NOTE: this does not check whether it's SE of PE--if you want to do that you'll have to handle it here.
        #tmp.append(input_template.format(sample=sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def getNormal_fastq(wildcards):
    return runsHelper(wildcards, 0)

def getTumor_fastq(wildcards):
    return runsHelper(wildcards, 1)



def somaticall_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #TODO- fill this in
    return ls

rule somaticcalls_all:
    input:
        somaticall_targets

rule somatic_calling_TNsnv:
    input:
        norm=getNormal,
        tumor=getTumor,
        corealignedbam="analysis/corealignments//{run}/{run}_tn_corealigned.bam"
    output:
        statscall="analysis/somaticVariants/{run}/{run}_call.output.stats"
        tnsnvvcf="analysis/somaticVariants/{run}/{run}_tnsnv.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads:_somaticcall_threads
    shell:
       ##min tumor allele fraction needs to be changed based on different value cutoffs##
       """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo TNsnv --tumor_sample {input.tumor} --normal_sample {input.normal} --dbsnp {params.dbsnp} --call_stats_out {output.statscall} --min_tumor_allele_frac 0.005 {output.tnsvvcf}"""


rule somatic_calling_TNhaplotyper:
    input:
        norm=getNormaL,
        tumor=getTumor,
        corealignedbam="analysis/corealignments//{run}/{run}_tn_corealigned.bam"
   output:
        tnhaplotypervcf="analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.vcf.gz"
   params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads:_somaticcall_threads
    shell:
       """{params.index1}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNhaplotyper --tumor_sample {input.tumor} --normal_sample {input.normal} --dbsnp {params.dbsnp} {output.tnhaplotypervcf}"""


rule somatic_calling_TNscope:
    input:
        norm=getNormaL,
        tumor=getTumor,
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
   output:
        tnscopevcf="analysis/somaticVariants/{run}/{run}_tnscope.output.vcf.gz"
   params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads:_somaticcall_threads
    shell:
       """{params.index1}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNscope --tumor_sample {input.tumor} --normal_sample {input.normal} --dbsnp {params.dbsnp} {output.tnscopevcf}"""

rule tnsnv_vcftoolsfilter:
    input:
        tnsnvvcf="analysis/somaticVariants/{run}/{run}_tnsnv.output.vcf.gz"
    output:
        tnsnvfilteredvcf="analysis/somaticVariants/{run}/{run}_tnsnv.output.filter.vcf"
    params:
        index=config['genome_fasta'],
    shell:
       """vcftools --gzvcf {input.tnsnvvcf} --remove-filtered-all --recode --stdout > {output.tnsnvfilteredvcf}"""

rule tnhaplotyper_vcftoolsfilter:
    input:
        tnhaplotypervcf="analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.vcf.gz"
    output:
        tnhaplotyperfilteredvcf="analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.filter.vcf"
    params:
        index=config['genome_fasta'],
    shell:
       """vcftools --gzvcf {input.tnhaplotypervcf} --remove-filtered-all --recode --stdout > {output.tnhaplotyperfilteredvcf}"""

rule tnscope_vcftoolsfilter:
    input:
        tnscopevcf="analysis/somaticVariants/{run}/{run}_tnscope.output.vcf.gz"
    output:
        tnscopefilteredvcf="analysis/somaticVariants/{run}/{run}_tnscope.output.filter.vcf"
    params:
        index=config['genome_fasta'],
    shell:
       """vcftools --gzvcf {input.tnscopevcf} --remove-filtered-all --recode --stdout > {output.tnscopefilteredvcf}"""
