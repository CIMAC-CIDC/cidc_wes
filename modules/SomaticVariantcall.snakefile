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
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.maf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.maf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnscope.output.maf" % (run,run))
        #EXON mutations
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.exon.maf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.exon.maf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnscope.output.exon.maf" % (run,run))

        #Mutation Signatures
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.pdf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.pdf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnscope.output.pdf" % (run,run))

        #alleleFrac cutoffs
        for frac in [0.05,0.1,0.2,0.3,0.4,0.5]:
            ls.append("analysis/somaticVariants/%s/%s_tnscope.output.%s.vcf" % (run,run, str(frac)))
            ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.%s.vcf" % (run,run, str(frac)))

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

rule vcftoolsfilter:
    """General rule to filter the three different types of vcf.gz files"""
    input:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.vcf.gz"
    output:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.filter.vcf"
    params:
        index=config['genome_fasta'],
    shell:
       """vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""


rule gunzip_vcf:
    """General rule to gunzip the three types of vcf.gz files-
    tnscope_, tnsnv, and tnhaplotyper"""
    input:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.vcf.gz"
    output:
        #Should we make this a temp?
        "analysis/somaticVariants/{run}/{run}_{caller}.output.vcf"
    shell:
        #NOTE: we want to keep the original .gz vcf file
        "gunzip -k {input}"
        
rule vcf2maf:
    """General rule to convert the different vcf files into maf"""
    input:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.vcf"
    output:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.maf"
    threads: _vcf2maf_threads,
    params:
        index=config['genome_fasta'],
        vep_path="%s/bin" % config['wes_root'],
        vep_data=config['vep_data'],
        vep_assembly=config['vep_assembly'],
    shell:
        """vcf2maf.pl --input-vcf {input} --output-maf {output} --ref-fasta {params.index} --vep-path {params.vep_path} --vep-data {params.vep_data} --ncbi-build {params.vep_assembly}"""  

rule mutationSignature:
    """General rule to do mutation signature analysis using mutProfiler.py"""
    input:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.maf"
    output:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.pdf"
    params:
        index= lambda wildcards: os.path.abspath(config['genome_fasta']),
        matrix="cidc_wes/cidc-vs/cidcvs/data/REF/TCGA-LUAD.mtrx", #HARD coding this for now!!!
        outname = lambda wildcards: "analysis/somaticVariants/%s/%s_%s.output" % (wildcards.run, wildcards.run, wildcards.caller),
        name = lambda wildcards: wildcards.run
    shell:
        "cidc_wes/cidc-vs/mutProfile.py -c {params.matrix} -m {input} -r {params.index} -o {params.outname} -n {params.name}"

rule alleleFrac_filter_tnscope:
    input:
        "analysis/somaticVariants/{run}/{run}_tnscope.output.vcf.gz"
    params:
        threshold=lambda wildcards: wildcards.frac
    output:
        #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
        #with vcftoolsfilter
        "analysis/somaticVariants/{run}/{run}_tnscope.output.{frac,\d\.\d+}.vcf"
    shell:
        "cidc_wes/modules/scripts/vcf_alleleFracFilter.py -v {input} -t {params.threshold} -o {output}"

rule alleleFrac_filter_tnhaplotyper:
    input:
        "analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.vcf.gz"
    params:
        threshold=lambda wildcards: wildcards.frac
    output:
        #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
        #with vcftoolsfilter
        "analysis/somaticVariants/{run}/{run}_tnhaplotyper.output.{frac,\d\.\d+}.vcf"
    shell:
        "cidc_wes/modules/scripts/vcf_alleleFracFilter.py -v {input} -t {params.threshold} -o {output}"

rule maf_exon_filter:
    """General rule to filter coding exon mutations"""
    input:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.maf"
    output:
        "analysis/somaticVariants/{run}/{run}_{caller}.output.exon.maf"
    shell:
        "cidc_wes/modules/scripts/maf_exon_filter.py -m {input} -o {output}"
