# module: Germline variant caller
# QC module to ensure that the tumor and normal samples come from the same
# patient

_germline_threads=32

def germline_runsHelper(wildcards, iindex, input_template):
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
    
def germline_getNormalTargets(wildcards):
    return germline_runsHelper(wildcards, 0, "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz")

def germline_getTumorTargets(wildcards):
    return germline_runsHelper(wildcards, 1, "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz")

def germline_getNormalTargetsTbi(wildcards):
    return germline_runsHelper(wildcards, 0, "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz.tbi")

def germline_getTumorTargetsTbi(wildcards):
    return germline_runsHelper(wildcards, 1, "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz.tbi")


def germline_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/germline/%s/%s_haplotyper.output.vcf" % (sample,sample))
        ls.append("analysis/germline/%s/%s_haplotyper.targets.vcf.gz" % (sample,sample))
        ls.append("analysis/germline/%s/%s_haplotyper.targets.vcf.gz.tbi" % (sample,sample))
    for run in config['runs']:
        ls.append("analysis/germline/%s/%s_vcfcompare.txt" % (run,run))
    return ls

def getTargetBed(config):
    """USES center_targets in somtatic.snakefile to return the path to the
    center's targets"""
    
    if 'cimac_center' in config and config['cimac_center'] in center_targets:
        center = config['cimac_center']
        return center_targets[center]
    else:
        return center_targets['broad']

rule germline_all:
    input:
        germline_targets

rule germline_haplotyper:
    input:
        recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam"
    output:
        haplotypervcf="analysis/germline/{sample}/{sample}_haplotyper.output.vcf"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
    threads:_germline_threads
    benchmark:
        "benchmarks/germline/{sample}/{sample}.germline_haplotyper.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.recalibratedbam} --algo Haplotyper  --dbsnp {params.dbsnp}  --emit_conf=30 --call_conf=30 {output.haplotypervcf}"""

rule germline_center_targets:
    """outputs the variants that are found only in the cimac_center's target
    bed file"""
    input:
        "analysis/germline/{sample}/{sample}_haplotyper.output.vcf"
    output:
        "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf"
    params:
        target_bed= lambda wildcards: getTargetBed(config['cimac_center'])
    benchmark:
        "benchmarks/germline/{sample}/{sample}.germline_center_targets.txt"    
    shell:
        """vcftools --vcf {input} --bed {params.target_bed} --recode --stdout > {output}"""

rule germline_bgzip:
    input:
        "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf"
    output:
        "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz"
    benchmark:
        "benchmarks/germline/{sample}/{sample}.germline_bgzip.txt"
    shell:
        "bgzip {input}"

rule germline_tabix:
    input:
        "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz"
    output:
        "analysis/germline/{sample}/{sample}_haplotyper.targets.vcf.gz.tbi"
    benchmark:
        "benchmarks/germline/{sample}/{sample}.germline_tabix.txt"
    shell:
        "tabix -p vcf {input}"
    
rule germline_vcfcompare:
    input:
        normal=germline_getNormalTargets,
        normal_tbi=germline_getNormalTargetsTbi,
        tumor=germline_getTumorTargets,
        turmor_tbi=germline_getTumorTargetsTbi,
    output:
        "analysis/germline/{run}/{run}_vcfcompare.txt",
    benchmark:
        "benchmarks/germline/{run}/{run}.germline_vcfcompare.txt"
    shell:
        "vcf-compare {input.tumor} {input.normal} > {output}"

#RULE to check {run}_vcfcompare.txt > 90% here
