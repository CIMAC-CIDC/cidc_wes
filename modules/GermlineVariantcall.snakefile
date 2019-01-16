#MODULE: Germline calls by Sentieon
#import os
#from string import Template

_germlinecall_threads=16

def germlinecalls_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/germlineVariants/%s/%s_dnascope.output.vcf.gz" % (run,run))
        ls.append("analysis/germlineVariants/%s/%s_haplotyper.output.vcf.gz" % (run,run))
    return ls

rule germlinecalls_all:
    input:
        germlinecalls_targets

rule germline_calling_DNAscope:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        dnascopevcf="analysis/germlineVariants/{run}/{run}_dnascope.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads:_germlinecall_threads
    shell:
       """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo DNAscope  --dbsnp {params.dbsnp}  --emit_conf=30 --call_conf=30 {output.dnascopevcf}"""

rule germline_calling_Haplotyper:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        haplotypervcf="analysis/germlineVariants/{run}/{run}_haplotyper.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads:_germlinecall_threads
    shell:
       """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo Haplotyper  --dbsnp {params.dbsnp} --emit_conf=30 --call_conf=30  {output.haplotypervcf}"""

