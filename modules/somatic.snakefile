#module: Somatic Variant calls by Sentieon
#import os
#from string import Template

_somatic_threads=32
#_vcf2maf_threads=4

#Dictionary of center targets
center_targets={'mocha':"./ref_files/hg38/target_beds/mocha.liftover.hg38.bed",
                "mda": "./ref_files/hg38/target_beds/MDA.liftover.hg38.bed",
                "broad":"./ref_files/hg38/target_beds/broad.liftover.hg38.bed"}

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


def somatic_getNormal(wildcards):
    return somatic_runsHelper(wildcards, 0)

def somatic_getTumor(wildcards):
    return somatic_runsHelper(wildcards, 1)

def somatic_getNormal_recal(wildcards):
    sample = somatic_runsHelper(wildcards, 0)[0]
    return "analysis/align/%s/%s_recalibrated.bam" % (sample,sample)

def somatic_getNTumor_recal(wildcards):
    sample = somatic_runsHelper(wildcards, 1)[0]
    return "analysis/align/%s/%s_recalibrated.bam" % (sample,sample)

def somatic_getNormal_recal_bai(wildcards):
    sample = somatic_runsHelper(wildcards, 0)[0]
    return "analysis/align/%s/%s_recalibrated.bam.bai" % (sample,sample)

def somatic_getNTumor_recal_bai(wildcards):
    sample = somatic_runsHelper(wildcards, 1)[0]
    return "analysis/align/%s/%s_recalibrated.bam.bai" % (sample,sample)

def somatic_tnsnv_targets(wildcards):
    ls = []
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s_tnsnv.output.vcf.gz" % (run,run))    
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.exons.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.output.vep.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.output.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.output.exon.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.output.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.mutationload.txt" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.stats.txt" % (run,run))
        for center in center_targets:
            ls.append("analysis/somatic/%s/%s_tnsnv.filter.exons.%s.vcf.gz" % (run,run,center))

    return ls

def somatic_tnhaplotyper2_targets(wildcards):
    ls = []
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.output.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.filter.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.filter.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.filter.exons.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.output.vep.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.output.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.output.exon.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.output.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.mutationload.txt" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper2.filter.stats.txt" % (run,run))
        for center in center_targets:
            ls.append("analysis/somatic/%s/%s_tnhaplotyper2.filter.exons.%s.vcf.gz" % (run,run,center))
    return ls

def somatic_tnscope_targets(wildcards):
    ls = []
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s_tnscope.output.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.exons.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.vep.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.exon.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.mutationload.txt" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.stats.txt" % (run,run))
        for center in center_targets:
            ls.append("analysis/somatic/%s/%s_tnscope.filter.exons.%s.vcf.gz" % (run,run,center))
    return ls

def somatic_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    if 'somatic_caller' in config:
        if config['somatic_caller'] == 'tnsnv':
            ls = somatic_tnsnv_targets(wildcards)
        elif config['somatic_caller'] == 'tnscope':
            ls = somatic_tnscope_targets(wildcards)
        else: #default
            ls = somatic_tnhaplotyper2_targets(wildcards)
    else: #default
        ls = somatic_tnhaplotyper2_targets(wildcards)
        
    #target summaries
    ls.append("analysis/metrics/all_sample_summaries.txt")
    return ls

rule somatic_all:
    input:
        somatic_targets

###############################################################################
# PLEASE LOOK AT the different caller snakefiles,
# somatic_tnhaplotyper2, somatic_tnsnv, somatic_tnscope
# for caller specific rules!
###############################################################################

rule filter_raw_vcf:
    """General rule to filter the three different types of vcf.gz files"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        tumor=lambda wildcards: somatic_getTumor(wildcards),
        normal= lambda wildcards: somatic_getNormal(wildcards),
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_filter_raw_vcf.txt"
    run: #DISABLES the conda env
        #SWITCH for tnhaplotyper2 filter
        if (wildcards.caller == "tnhaplotyper2"):
            shell("{params.sentieon_path}/sentieon tnhapfilter --tumor_sample {params.tumor} --normal_sample {params.normal} -v {input} {output}")
        else:
            shell("""vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}""")


rule gunzip_vcf:
    """General rule to gunzip the three types of vcf.gz files-
    tnscope_, tnsnv, and tnhaplotyper"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz"
    output:
        #Should we make this a temp?
        "analysis/somatic/{run}/{run}_{caller}.output.vcf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_gunzip_vcf.txt"
    group: "somatic"
    shell:
        #NOTE: we want to keep the original .gz vcf file
        "gunzip < {input} > {output}"

rule vcfVEP:
    """Rule to annotate vcf files with vep"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.{type}.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.vep.vcf"
    params:
        vep_data=config['vep_data'],
        vep_synonyms=config['vep_synonyms'],
        gdc_fasta=config['genome_fasta'],
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcfVEP.txt"
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        "vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs --fa {params.gdc_fasta}"
    
rule vcf2maf:
    """General rule to convert the different vcf files into maf"""
    input:
        vcf="analysis/somatic/{run}/{run}_{caller}.{type}.vcf",
        vep="analysis/somatic/{run}/{run}_{caller}.{type}.vep.vcf",
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.maf"
    params:
        vep_index=config['vep_fasta'],
        vep_custom_enst= config['vep_custom_enst'],
        vep_assembly=config['vep_assembly'],
        vep_filter= config['vep_filter'],
        buffer_size=config['vcf2maf_bufferSize'],

        tumor= lambda wildcards: somatic_getTumor(wildcards),
        normal= lambda wildcards: somatic_getNormal(wildcards),
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcf2maf.txt"
    log:
        "analysis/logs/somatic/{run}/{run}.{caller}.{type}_vcf2maf.log.txt"
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        """vcf2maf.pl --input-vcf {input} --output-maf {output} --custom-enst {params.vep_custom_enst} --ref-fasta {params.vep_index} --tumor-id {params.tumor} --normal-id {params.normal} --ncbi-build {params.vep_assembly} --filter-vcf {params.vep_filter} --buffer-size {params.buffer_size} 2> {log}"""


rule mutationSignature:
    """General rule to do mutation signature analysis using mutProfiler.py"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.{type}.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.pdf"
    params:
        index= lambda wildcards: os.path.abspath(config['genome_fasta']),
        matrix="cidc_wes/cidc-vs/cidcvs/data/REF/TCGA-LUAD.mtrx", #HARD coding this for now!!!
        outname = lambda wildcards: "%sanalysis/somatic/%s/%s_%s.%s" % (config['remote_path'], wildcards.run, wildcards.run, wildcards.caller, wildcards.type),
        name = lambda wildcards: wildcards.run
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_mutationSignature.txt"
    group: "somatic"
    shell:
        "cidc_wes/cidc-vs/mutProfile.py -c {params.matrix} -m {input} -r {params.index} -o {params.outname} -n {params.name}"

rule maf_exon_filter:
    """General rule to filter coding exon mutations"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_maf_exon_filter.txt"
    shell:
        "cidc_wes/modules/scripts/maf_exon_filter.py -m {input} -o {output}"

rule calculate_mutation:
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.mutationload.txt"
    params:
        size=config['effective_size'],
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.calculate_mutation.txt"
    shell:
        "cidc_wes/modules/scripts/mutation_load.py -v {input} -o {output} -s {params.size}"

rule extract_VAF_DEPTH:
    """Run Jingxins harmonization script to extract VAF and DEPTH stats"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.stats.txt"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.extract_VAF_DEPTH.txt"
    shell:
        "cidc_wes/modules/scripts/extract_vaf_depth.py -v {input} > {output}"

rule somatic_gzip_filtered_vcf:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz",
    group: "somatic"
    shell:
        "bgzip -c {input} > {output}"

rule somatic_tabix_filtered_vcf_gz:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz.tbi",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        "tabix -p vcf {input}"

rule somatic_getExonic_mutations:
    """Get the mutations that fall into the exonic regions"""
    input:
        vcf="analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz",
        tbi="analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz.tbi"
    params:
        exons=config['CDS_Bed_input']
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.somatic_getExonic_mutations.txt"
    shell:
        "bcftools view -R {params.exons} {input.vcf} | bcftools sort | bcftools view -Oz > {output}"

rule somatic_tabix_exonic_mutations:
    """Get the mutations that fall into the exonic regions"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz.tbi",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        "tabix -p vcf {input}"

rule somatic_getTarget_mutations:
    """Get the mutations that fall into the exonic regions"""
    input:
        vcf="analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
        tbi="analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz.tbi"
    params:
        target= lambda wildcards: center_targets[wildcards.center]
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.{center}.vcf.gz",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.{center}.somatic_getTarget_mutations.txt"
    shell:
        "bcftools view -R {params.target} {input} | bcftools sort | bcftools view -Oz > {output}"

rule somatic_collect_target_summaries:
    """Collect all of the sample summaries and put them in one file
    input: {sample}_target_metrics.txt.sample_summary from coverage.snakefile
    """
    input:
        expand("analysis/metrics/{sample}/{sample}_target_metrics.txt.sample_summary", sample=sorted(config['samples']))
    output:
        "analysis/metrics/all_sample_summaries.txt"
    params:
        files = lambda wildcards, input: " -f ".join(input)
    group: "somatic"
    benchmark:
        "benchmarks/somatic/somatic_collect_target_summaries.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_collect_target_summaries.py -f {params.files} > {output}"
