#module: Somatic Variant calls by Sentieon
#import os
#from string import Template

_somaticcall_threads=32
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

def somaticall_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #Consolidate these with an inner-for-loop?
        ls.append("analysis/somatic/%s/%s_call.output.stats" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnsnv.output.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.output.vcf.gz" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.vcf.gz" % (run,run))
        #FILTERED VCF
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.filter.vcf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.vcf" % (run,run))
        #MAF
        ls.append("analysis/somatic/%s/%s_tnsnv.output.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.output.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.maf" % (run,run))
        #Filtered MAF
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.filter.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.maf" % (run,run))
        #Mutation Signatures
        ls.append("analysis/somatic/%s/%s_tnsnv.output.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.output.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.pdf" % (run,run))
        #Filtered Mutation Signatures
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.filter.pdf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.filter.pdf" % (run,run))

        #EXON mutations- should this be on full or filtered?
        ls.append("analysis/somatic/%s/%s_tnsnv.output.exon.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnhaplotyper.output.exon.maf" % (run,run))
        ls.append("analysis/somatic/%s/%s_tnscope.output.exon.maf" % (run,run))
        #alleleFrac cutoffs - should this be on full or filtered?
        for frac in [0.05,0.1,0.2,0.3,0.4,0.5]:
            ls.append("analysis/somatic/%s/%s_tnscope.output.%s.vcf" % (run,run, str(frac)))
            ls.append("analysis/somatic/%s/%s_tnhaplotyper.output.%s.vcf" % (run,run, str(frac)))

        #read depth/coverage filter: 10x, 20x, 50x - should this be on full or filtered?
        for frac in [10, 20, 50]:
            ls.append("analysis/somatic/%s/%s_tnscope.coverage.%s.vcf" % (run,run, str(frac)))
            ls.append("analysis/somatic/%s/%s_tnsnv.coverage.%s.vcf" % (run,run, str(frac)))

        #BEGIN HARMONIZATION
        #Mutation load
        ls.append("analysis/somatic/%s/%s_tnsnv.mutationload.txt" % (run,run))
        #STATS:
        ls.append("analysis/somatic/%s/%s_tnsnv.filter.stats.txt" % (run,run))
        #Center specific exon targets
        for center in center_targets:
            ls.append("analysis/somatic/%s/%s_tnsnv.filter.exons.%s.vcf.gz" % (run,run,center))
    return ls

rule somaticcalls_all:
    input:
        somaticall_targets

rule somatic_calling_TNsnv:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        statscall="analysis/somatic/{run}/{run}_call.output.stats",
        tnsnvvcf="analysis/somatic/{run}/{run}_tnsnv.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        #JUST sample names - can also use the helper fns, e.g.
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:96
    benchmark:
        "benchmarks/somatic/{run}/{run}.somatic_calling_TNsnv.txt"
    shell:
        #"""{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo TNsnv --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} --call_stats_out {output.statscall} --min_tumor_allele_frac 0.05 {output.tnsnvvcf}"""
        #REMOVING min_tumor_allele_frac param
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo TNsnv --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} --call_stats_out {output.statscall} {output.tnsnvvcf}"""


rule somatic_calling_TNhaplotyper:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        tnhaplotypervcf="analysis/somatic/{run}/{run}_tnhaplotyper.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: getNormal_sample(wildcards)
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:_somaticcall_threads
    benchmark:
        "benchmarks/somatic/{run}/{run}.somatic_calling_TNhaplotyper.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNhaplotyper --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {output.tnhaplotypervcf}"""


rule somatic_calling_TNscope:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        tnscopevcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        #JUST sample names - can also use the helper fns, e.g.
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:_somaticcall_threads
    benchmark:
        "benchmarks/somatic/{run}/{run}.somatic_calling_TNscope.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNscope --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {output.tnscopevcf}"""

rule vcftoolsfilter:
    """General rule to filter the three different types of vcf.gz files"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    params:
        index=config['genome_fasta'],
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_vcftoolsfilter.txt"
    shell:
       """vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""

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
    shell:
        #NOTE: we want to keep the original .gz vcf file
        "gunzip -k {input}"

rule vcfVEP:
    """Rule to annotate vcf files with vep"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.{type}.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.vep.vcf"
    params:
        vep_data=config['vep_data'],
        vep_synonyms=config['vep_synonyms'],
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcfVEP.txt"
    shell:
        "vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs"
    
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
        "analysis/log/somaticvariantcall/{run}/{run}.{caller}.{type}_vcf2maf.log.txt"
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
        outname = lambda wildcards: "analysis/somatic/%s/%s_%s.%s" % (wildcards.run, wildcards.run, wildcards.caller, wildcards.type),
        name = lambda wildcards: wildcards.run
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_mutationSignature.txt"
    shell:
        "cidc_wes/cidc-vs/mutProfile.py -c {params.matrix} -m {input} -r {params.index} -o {params.outname} -n {params.name}"

rule alleleFrac_filter_tnscope:
    input:
        "analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
    params:
        threshold=lambda wildcards: wildcards.frac
    output:
        #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
        #with vcftoolsfilter
        "analysis/somatic/{run}/{run}_tnscope.output.{frac,\d\.\d+}.vcf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.alleleFrac_filter_tnscope.txt"
    shell:
        "cidc_wes/modules/scripts/vcf_alleleFracFilter.py -v {input} -t {params.threshold} -o {output}"

rule alleleFrac_filter_tnhaplotyper:
    input:
        "analysis/somatic/{run}/{run}_tnhaplotyper.output.vcf.gz"
    params:
        threshold=lambda wildcards: wildcards.frac
    output:
        #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
        #with vcftoolsfilter
        "analysis/somatic/{run}/{run}_tnhaplotyper.output.{frac,\d\.\d+}.vcf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.alleleFrac_filter_tnhaplotyper.txt"
    shell:
        "cidc_wes/modules/scripts/vcf_alleleFracFilter.py -v {input} -t {params.threshold} -o {output}"

rule maf_exon_filter:
    """General rule to filter coding exon mutations"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_maf_exon_filter.txt"
    shell:
        "cidc_wes/modules/scripts/maf_exon_filter.py -m {input} -o {output}"

rule coverage_filter_tnscope:
    input:
        "analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
    params:
        threshold=lambda wildcards: wildcards.frac,
        field="AFDP" #NOTE this is the particular field tnscope vcf files
    output:
        #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
        #with vcftoolsfilter; {frac} is int
        "analysis/somatic/{run}/{run}_tnscope.coverage.{frac,\d+}.vcf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.coverage_filter_tnscope.txt"
    shell:
        "cidc_wes/modules/scripts/vcf_filterByReadDepth.py -v {input} -t {params.threshold} -f {params.field} -o {output}"

rule coverage_filter_tnsnv:
    input:
        "analysis/somatic/{run}/{run}_tnsnv.output.vcf.gz"
    params:
        threshold=lambda wildcards: wildcards.frac,
        field="DP" #NOTE this is the particular field tnsnv vcf files
    output:
        #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
        #with vcftoolsfilter; {frac} is int
        "analysis/somatic/{run}/{run}_tnsnv.coverage.{frac,\d+}.vcf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.coverage_filter_tnsnv.txt"
    shell:
        "cidc_wes/modules/scripts/vcf_filterByReadDepth.py -v {input} -t {params.threshold} -f {params.field} -o {output}"

rule calculate_mutation:
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.mutationload.txt"
    params:
        size=config['effective_size'],
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
    shell:
        "bgzip -c {input} > {output}"

rule somatic_tabix_filtered_vcf_gz:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz.tbi",
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
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.somatic_getTarget_mutations.txt"
    shell:
        "bcftools view -R {params.target} {input} | bcftools sort | bcftools view -Oz > {output}"