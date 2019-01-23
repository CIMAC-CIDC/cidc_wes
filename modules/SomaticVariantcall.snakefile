#MODULE: Somatic Variant calls by Sentieon
#import os
#from string import Template

_somaticcall_threads=16
_vcf2maf_threads=4

#NOTE: somatic_runsHelper, getNormal_sample, and getTumor_sample are NOT
#called by any one!
def somatic_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append(sample_name) 
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def getNormal_sample(wildcards):
    return somatic_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return somatic_runsHelper(wildcards, 1)

def somaticall_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/somaticVariants/%s/%s_call.output.stats" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.vcf.gz" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.vcf.gz" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnscope.output.vcf.gz" % (run,run))
        #FILTERED VCF
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.filter.vcf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.filter.vcf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnscope.output.filter.vcf" % (run,run))
        #MAF
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.filter.maf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.filter.maf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnscope.output.filter.maf" % (run,run))
    return ls

rule somaticcalls_all:
    input:
        somaticall_targets

rule somatic_calling_TNsnv:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        statscall="analysis/somaticVariants/{run}/{run}_call.output.stats",
        tnsnvvcf="analysis/somaticVariants/{run}/{run}_tnsnv.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: getNormal_sample(wildcards)
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:_somaticcall_threads
    shell:
        ##min tumor allele fraction needs to be changed based on different value cutoffs##
        #LT: will need to iterate over this list of values and run shell each time
        ##Values={0.05,0.1,0.2,0.3,0.4,0.5}
        #LT: will need to iterate over this list of values and run shell each time
        """{params.index1}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo TNsnv --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} --call_stats_out {output.statscall} --min_tumor_allele_frac 0.05 {output.tnsnvvcf}"""


rule somatic_calling_TNhaplotyper:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        tnhaplotypervcf="analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: getNormal_sample(wildcards)
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:_somaticcall_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNhaplotyper --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {output.tnhaplotypervcf}"""


rule somatic_calling_TNscope:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        tnscopevcf="analysis/somaticVariants/{run}/{run}_tnscope.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: getNormal_sample(wildcards)
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:_somaticcall_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNscope --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {output.tnscopevcf}"""

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

rule tnsnv_vcf2maf:
    input:
        tnsnvvcf="analysis/somaticVariants/{run}/{run}_tnsnv.output.filter.vcf"
    output:
        tnsnvmaf="analysis/somaticVariants/{run}/{run}_tnsnv.output.filter.maf"
    threads: _vcf2maf_threads,
    params:
        index=config['genome_fasta'],
        vep_path="%s/bin" % config['wes_root'],
        vep_data=config['vep_data'],
        vep_assembly=config['vep_assembly'],
    shell:
        #""" zcat {input.tnsnvvcf}; perl vcf2maf.pl --input-vcf - --output-maf {output.tnsnvmaf} --ref-fasta {params.index}"""
        """vcf2maf.pl --input-vcf {input.tnsnvvcf} --output-maf {output.tnsnvmaf} --ref-fasta {params.index} --vep-path {params.vep_path} --vep-data {params.vep_data} --ncbi-build {params.vep_assembly}"""  

rule tnhaplotyper_vcf2maf:
    input:
        tnhaplotypervcf="analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.filter.vcf"
    output:
        tnhaplotypermaf="analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.filter.maf"
    threads: _vcf2maf_threads,
    params:
        index=config['genome_fasta'],
        vep_path="%s/bin" % config['wes_root'],
        vep_data=config['vep_data'],
        vep_assembly=config['vep_assembly'],
    shell:
        #""" zcat {input.tnhaplotypervcf}; perl vcf2maf.pl --input-vcf - --output-maf {output.tnhaplotypermaf} --ref-fasta {params.index}"""
        """vcf2maf.pl --input-vcf {input.tnhaplotypervcf} --output-maf {output.tnhaplotypermaf} --ref-fasta {params.index}  --vep-path {params.vep_path} --vep-data {params.vep_data} --ncbi-build {params.vep_assembly}"""

rule tnscope_vcf2maf:
    input:
        tnscopevcf="analysis/somaticVariants/{run}/{run}_tnscope.output.filter.vcf"
    output:
        tnscopemaf="analysis/somaticVariants/{run}/{run}_tnscope.output.filter.maf"
    threads: _vcf2maf_threads,
    params:
        index=config['genome_fasta'],
        vep_path="%s/bin" % config['wes_root'],
        vep_data=config['vep_data'],
        vep_assembly=config['vep_assembly'],
    shell:
        #""" zcat {input.tnscopevcf}; perl vcf2maf.pl --input-vcf - --output-maf {output.tnscopemaf} --ref-fasta {params.index}"""
        """vcf2maf.pl --input-vcf {input.tnscopevcf} --output-maf {output.tnscopemaf} --ref-fasta {params.index}  --vep-path {params.vep_path} --vep-data {params.vep_data} --ncbi-build {params.vep_assembly}"""


