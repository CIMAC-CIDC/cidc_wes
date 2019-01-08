#MODULE: Germline calls by Sentieon
#import os
#from string import Template

_germlinecall_threads=8
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

def germlinecalls_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #TODO- fill this in
    return ls

rule germlinecalls_all:
    input:
        germlinecalls_targets

rule germline_calling_DNAscope:
    input:
        norm=getNormal,
        realignedbam="analysis/align/{sample}/{sample}_recalibrated.bam"
    output:
        dnascopevcf="analysis/germlineVariants/{run}/{run}_dnascope.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads:_germlinecall_threads
    shell:
       """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo DNAscope  --dbsnp {params.dbsnp}  --emit_conf=30 --call_conf=30 {output.dnascopevcf}"""

rule germline_calling_Haplotyper:
    input:
        norm=getNormal,
        realignedbam="analysis/align/{sample}/{sample}_recalibrated.bam"
    output:
        haplotypervcf="analysis/germlineVariants/{run}/{run}_haplotyper.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp=_dbsnp,
        mills=_Mills_indels,
        g1000=_G1000_indels,
    threads:_germlinecall_threads
    shell:
       """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo Haplotyper  --dbsnp {params.dbsnp} --emit_conf=30 --call_conf=30  {output.haplotypervcf}"""

