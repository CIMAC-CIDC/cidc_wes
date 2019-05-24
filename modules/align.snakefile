#MODULE: Align fastq files to genome - common rules
#import os
#_logfile="analysis/logs/align.log"
_align_threads=32
_bwa_threads=16

def align_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.bam.bai" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.dedup.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.dedup.bam.bai" % (sample,sample))
    ls.append("analysis/align/mapping.csv")
    return ls

def align_mapping_targets(wildcards):
    """Generates just the mapping targets for rule map_all"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (sample,sample))
    return ls

def align_getFastq(wildcards):
    ls = config["samples"][wildcards.sample]
    return ls

rule align_all:
    input:
        align_targets

rule map_all:
    input:
        align_mapping_targets

rule sentieon_bwa:
    input:
        align_getFastq
    output:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai"
    params:
        sentieon_path=config['sentieon_path'],
        bwa_index=config['bwa_index'],
        read_group= lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample),        
        input_bases="10000000",
        #need to adjust threads for the other process
        tthreads=lambda wildcards, input, output, threads, resources: threads-1
    threads: _bwa_threads
    message: "ALIGN: Running sentieon BWA mem for alignment"
    log: "analysis/logs/align.sentieon_bwa.{sample}.log"
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.sentieon_bwa.txt"
    shell:
        """({params.sentieon_path}/sentieon bwa mem -M -R \"{params.read_group}\" -t {params.tthreads} -K {params.input_bases} {params.bwa_index} {input} || echo -n 'error' ) | {params.sentieon_path}/sentieon util sort -r {params.bwa_index} -o {output.bam} --sam2bam -i -"""

rule map_stats:
    """Get the mapping stats for each aligment run"""
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
    output:
        "analysis/align/{sample}/{sample}_mapping.txt"
    threads: _align_threads
    message: "ALIGN: get mapping stats for each bam"
    log: "analysis/logs/align.map_stats.{sample}.log"
    #CAN/should samtools view be multi-threaded--
    #UPDATE: tricky on how to do this right w/ compounded commands
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.map_stats.txt"
    shell:
        #FLAGSTATS is the top of the file, and we append the uniquely mapped
        #reads to the end of the file
        "sambamba flagstat -t {threads} {input.bam} > {output} 2>>{log}"
        #" && sambamba view -c -t {threads} {input.uniq_bam} >> {output} 2>> {log}"

rule collect_map_stats:
    """Collect and parse out the mapping stats for the ALL of the samples"""
    input:
        #samples sorted to match order of rest of report
        expand("analysis/align/{sample}/{sample}_mapping.txt", sample=sorted(config["samples"]))
    output:
        "analysis/align/mapping.csv"
    message: "ALIGN: collect and parse ALL mapping stats"
    log: "analysis/logs/align.collect_map_stats.log"
    group: "align"
    benchmark:
        "benchmarks/align/collect_map_stats.txt"
    run:
        files = " -f ".join(input)
        shell("cidc_wes/modules/scripts/align_getMapStats.py -f {files} > {output} 2>>{log}")

rule scoreSample:
    "Calls sentieon driver  --fun score_info on the sample"
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    output:
        temp("analysis/align/{sample}/{sample}.sorted.score.txt")
    message: "ALIGN: score sample"
    log: "analysis/logs/align.scoreSample.{sample}.log"
    threads: _align_threads
    params:
        index1=config['sentieon_path'],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.scoreSample.txt"
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output}"""

rule dedupSortedUniqueBam:
    """Dedup sorted unique bams using sentieon
     output {sample}_unique.sorted.dedup.bam"""
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
        score="analysis/align/{sample}/{sample}.sorted.score.txt"
    output:
        bamm="analysis/align/{sample}/{sample}.sorted.dedup.bam",
        baii="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
        met="analysis/align/{sample}/{sample}.sorted.dedup.metric.txt",
    message: "ALIGN: dedup sorted unique bam file"
    log: "analysis/logs/align.dedupSortedUniqueBam.{sample}.log"
    threads: _align_threads
    params:
        index1=config['sentieon_path'],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.dedupSortedUniqueBam.txt"
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {input.score} --metrics {output.met} {output.bamm}"""
