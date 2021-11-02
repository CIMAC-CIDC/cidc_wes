# author: Len Taing (TGBTG) 
# year: 2020
# module: Sentieon RNAseq variant calling pipeline
# ref: https://support.sentieon.com/manual/RNA_call/rna/

_rna_threads=8
#map from tumor samples to runs
_tumor_run_map = dict([(config['runs'][r][1], r) for r in config['runs']])

def rna_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    rna = config.get("rna", None)
    caller = config.get('somatic_caller', 'tnscope')
    if rna:
        for tumor_sample in rna:
            run = _tumor_run_map.get(tumor_sample, None)            
            if run:
                ls.append("analysis/rna/%s/%s.RG.dedup.split.bam" % (run, tumor_sample))
                ls.append("analysis/rna/%s/%s.RG.dedup.recal.csv" % (run, tumor_sample))
                ls.append("analysis/rna/%s/%s.haplotyper.rna.vcf.gz" % (run, run))
                ls.append("analysis/rna/%s/%s_%s.output.twist.neoantigen.vep.rna.vcf" % (run, run, caller))
    #print(ls)
    return ls

rule rna_all:
    input:
        rna_targets
    #benchmark: "benchmarks/rna/rna_all.txt"

def rna_addReadGroup_inputFn(wildcards):
    sample = config['rna'].get(wildcards.sample, None)
    return sample['bam_file']
    
rule rna_addReadGroup:
    """NOTE: RIMA's STAR files don't have read group info--inject the 
    wes sample info
    """
    input:
        rna_addReadGroup_inputFn
    output:
        "analysis/rna/{run}/{sample}.RG.bam"
    params:
        iid = lambda wildcards: "ID:%s" % wildcards.sample,
        sm = lambda wildcards: "SM:%s" % wildcards.sample,
        pl = "PL:ILLUMINA",
    threads: 4
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{sample}.rna_addReadGroup.txt"
    shell:
        """samtools addreplacerg -r {params.iid} -r {params.sm} -r {params.pl} -@ {threads} -o {output} {input}"""

rule rna_indexRGbam:
    input:
        "analysis/rna/{run}/{sample}.RG.bam"
    output:
        "analysis/rna/{run}/{sample}.RG.bam.bai"
    threads: 4 #_rna_threads
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{sample}.rna_indexRGbam.txt"
    shell:
        """sambamba index -t {threads} {input}"""


        
#------------------------------------------------------------------------------
# DEDUP rules
rule rna_scoreBam:
    input:
        bam="analysis/rna/{run}/{sample}.RG.bam",
        bai="analysis/rna/{run}/{sample}.RG.bam.bai",
    output:
        "analysis/rna/{run}/{sample}.RG.score.txt"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
    threads: 4 #_rna_threads
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{sample}.rna_scoreBam.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.bam} --algo LocusCollector --fun score_info {output}"""

rule rna_dedup:
    input:
        bam="analysis/rna/{run}/{sample}.RG.bam",
        bai="analysis/rna/{run}/{sample}.RG.bam.bai",
        score="analysis/rna/{run}/{sample}.RG.score.txt",
    output:
        bam="analysis/rna/{run}/{sample}.RG.dedup.bam",
        met="analysis/rna/{run}/{sample}.RG.dedup.metrics.txt"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
    threads: 6 #_rna_threads
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{sample}.rna_dedup.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.bam} --algo Dedup --rmdup --score_info {input.score} --metrics {output.met} {output.bam}"""
#------------------------------------------------------------------------------

rule rna_BQSR:
    input:
        "analysis/rna/{run}/{sample}.RG.dedup.bam"
    output:
        "analysis/rna/{run}/{sample}.RG.dedup.recal.csv"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        mills= config['Mills_indels'],
        g1000= config['G1000_indels'],
    threads: 4 #_rna_threads
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{sample}.rna_BQSR.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input} --algo QualCal -k {params.dbsnp} -k {params.mills} -k {params.g1000} {output}"""

rule rna_splitReadsAtJunct:
    input:
        "analysis/rna/{run}/{sample}.RG.dedup.bam"
    output:
        "analysis/rna/{run}/{sample}.RG.dedup.split.bam"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
    threads: 4 #_rna_threads
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{sample}.rna_splitReadsAtJunct.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {output}"""

def rna_variantCalling_inputFn(wildcards):
    run = wildcards.run
    sample = config['runs'][run][1]
    caller = config.get('somatic_caller', 'tnscope')
    tmp = {}
    tmp['bam'] = "analysis/rna/%s/%s.RG.dedup.split.bam" % (run, sample)
    tmp['bqsr'] = "analysis/rna/%s/%s.RG.dedup.recal.csv" % (run, sample)
    #USE this instead
    tmp['vcf'] = "analysis/somatic/%s/%s_%s.output.twist.neoantigen.vep.vcf.gz" % (run, run, caller)
    return tmp

rule rna_haplotyper:
    """Performs the sentieon haplotyper variant caller using the 
    filter.neoantigen.vep.vcf.gz file as a list of given variants to check.
    The shell call is the same as the rna pipeline found here:
    https://support.sentieon.com/manual/RNA_call/rna/
    """
    input:
        unpack(rna_variantCalling_inputFn)
    output:
        "analysis/rna/{run}/{run}.haplotyper.rna.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
    threads: 4 #_rna_threads
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{run}.rna_haplotyper.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} -q {input.bqsr}  --algo Haplotyper  --trim_soft_clip --call_conf 20 --emit_conf 20 -d {params.dbsnp} --given {input.vcf} {output}"""
        #DNASCOPE-
        #"""{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.bam} -q {input.bqsr}  --algo DNAscope --trim_soft_clip --call_conf 20 --emit_conf 20 -d {params.dbsnp} --given {input.vcf} {output}"""

def rna_intersect_inputFn(wildcards):
    run = wildcards.run
    caller = config.get('somatic_caller', 'tnscope')
    tmp = {}
    tmp['dna'] = "analysis/somatic/%s/%s_%s.output.twist.neoantigen.vep.vcf.gz" % (run, run, caller)
    tmp['rna'] = "analysis/rna/%s/%s.haplotyper.rna.vcf.gz" % (run, run)
    return tmp

rule rna_intersect:
    """Intersect WES somatic variants with RNA variants; 
    return record from WES"""
    input:
        unpack(rna_intersect_inputFn)
    output:
        "analysis/rna/{run}/{run}_{caller}.output.twist.neoantigen.vep.rna.vcf"
    group: "rna"
    benchmark:
        "benchmarks/rna/{run}/{run}_{caller}.rna_intersect.txt"
    shell:
        """bcftools isec {input.dna} {input.rna} -n =2 -w 1 -O v > {output}"""


    
