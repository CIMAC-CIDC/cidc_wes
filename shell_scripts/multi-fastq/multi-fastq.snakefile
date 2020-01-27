configfile: "config.yaml"
_bwa_threads=32

def targets(wildcards):
    ls = []
    for sample in config['samples']:
        #FOR MERGED
        #ls.append("analysis/merged/%s/%s.sorted.bam" % (sample, sample))

        #FOR Sentieon
        ls.append("analysis/dedup/%s/%s.sorted.dedup.bam" % (sample, sample))
    return ls

def align_getFastq(wildcards):
    ls = config["read_groups"][wildcards.read_group]
    return ls

def getSample(read_group):
    """Given a read group, returns the sample the RG belongs to (string)"""
    samples = config['samples']
    for s in samples:
        if read_group in config['samples'][s]:
            return s
    #CASE- read group not found, return ""
    return ""

def getReadGroups(wildcards):
    """Given a wildcards.sample, returns list of (aligned) read group files"""
    ls = []
    sample = wildcards.sample
    for rg in config['samples'][sample]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (rg,rg))
    return ls

rule all:
    input:
        targets

rule align_from_fastq:
    input:
        align_getFastq
    output:
        bam="analysis/align/{read_group}/{read_group}.sorted.bam",
        bai="analysis/align/{read_group}/{read_group}.sorted.bam.bai"
    params:
        sentieon_path=config['sentieon_path'],
        bwa_index=config['bwa_index'],
        #GET sample that the read_group belongs to using getSample
        RG= lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.read_group, getSample(wildcards.read_group)),
        input_bases="10000000",
    threads: _bwa_threads
    message: "ALIGN: Running sentieon BWA mem for alignment"
    log: "analysis/logs/align/{read_group}/align.sentieon_bwa.{read_group}.log"
    group: "align"
    benchmark:
        "benchmarks/align/{read_group}/{read_group}.align_from_fastq.txt"
    shell:
        """({params.sentieon_path}/sentieon bwa mem -M -R \"{params.RG}\" -t {threads} -K {params.input_bases} {params.bwa_index} {input} || echo -n 'error' ) | {params.sentieon_path}/sentieon util sort -r {params.bwa_index} -o {output.bam} --sam2bam -i -"""
    
rule merge:
    input:
        getReadGroups
    output:
        bam="analysis/merge/{sample}/{sample}.sorted.bam",
    params:
        sentieon_path=config['sentieon_path'],
    threads: _bwa_threads
    message: "MERGE: merging all sample parts"
    log: "analysis/logs/merge/{sample}/align.merge.{sample}.log"
    group: "align"
    benchmark:
        "benchmarks/merge/{sample}/{sample}.merge.txt"
    shell:
        "{params.sentieon_path}/sentieon util sort -t {threads} -o {output} {input}"
    
###############################################################################
# SENTIEON multi-fastq method below
###############################################################################
rule dedup_score_sample:
    input:
        getReadGroups
    output:
        score="analysis/dedup/{sample}/{sample}.sorted.score.txt",
        idx="analysis/dedup/{sample}/{sample}.sorted.score.txt.idx",
    message: "ALIGN: score sample"
    log: "analysis/logs/dedup/{sample}/align.scoreSample.{sample}.log"
    threads: _bwa_threads
    params:
        sentieon_path=config['sentieon_path'],
        inputs = lambda wildcards,input: [" -i %s" % i for i in input],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.scoreSample.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} {params.inputs} --algo LocusCollector --fun score_info {output.score}"""

rule dedup_bams:
    input:
        rg=getReadGroups,
        score_input="analysis/dedup/{sample}/{sample}.sorted.score.txt",
    output:
        bamm="analysis/dedup/{sample}/{sample}.sorted.dedup.bam",
        baii="analysis/dedup/{sample}/{sample}.sorted.dedup.bam.bai",
        met="analysis/dedup/{sample}/{sample}.sorted.dedup.metric.txt",
    message: "ALIGN: dedup sorted unique bam file"
    log: "analysis/logs/dedup/{sample}/align.dedupSortedUniqueBam.{sample}.log"
    threads: _bwa_threads
    params:
        sentieon_path=config['sentieon_path'],
        inputs = lambda wildcards,input: [" -i %s" % i for i in input.rg],

    group: "align"
    benchmark:
        "benchmarks/dedup/{sample}/{sample}.dedupSortedUniqueBam.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} {params.inputs} --algo Dedup --rmdup --score_info {input.score_input} --metrics {output.met} {output.bamm}"""

