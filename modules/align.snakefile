#MODULE: Align fastq files to genome - common rules
#import os
_logfile="analysis/logs/align.log"
_align_threads=32
_bwa_threads=16

def align_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.bam.bai" % (sample,sample))
        
        #REMOVING THIS!
        #ls.append("analysis/align/%s/%s_unique.sorted.bam" % (sample,sample))
        #ls.append("analysis/align/%s/%s_unique.sorted.bam.bai"%(sample,sample))
        ls.append("analysis/align/%s/%s.sorted.dedup.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.dedup.bam.bai" % (sample,sample))
    ls.append("analysis/align/mapping.csv")
    return ls

def all_mapping_targets(wildcards):
    """Generates just the mapping targets for rule map_all"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (sample,sample))
    return ls

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule align_all:
    input:
        align_targets

rule map_all:
    input:
        all_mapping_targets

rule sentieon_bwa:
    input:
        getFastq
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
    log: _logfile
    shell:
        """({params.sentieon_path}/sentieon bwa mem -M -R \"{params.read_group}\" -t {params.tthreads} -K {params.input_bases} {params.bwa_index} {input} || echo -n 'error' ) | {params.sentieon_path}/sentieon util sort -r {params.bwa_index} -o {output.bam} --sam2bam -i -"""

#REMOVING THIS
# rule uniquely_mapped_reads:
#     """Get the uniquely mapped reads"""
#     input:
#         "analysis/align/{sample}/{sample}.sorted.bam"
#     params:
#         filter="\'mapping_quality >= 1\'"
#     output:
#         "analysis/align/{sample}/{sample}_unique.sorted.bam"
#     message: "ALIGN: Filtering for uniquely mapped reads"
#     log: _logfile
#     threads: _align_threads
#     shell:
#         #NOTE: this is the generally accepted way of doing this as multiply 
#         #mapped reads have a Quality score of 0
#         #"samtools view -bq 1 -@ {threads} {input} > {output}"
#         "sambamba view -f bam -F {params.filter} -t {threads} {input} > {output}"

#NOTE: dropping uniquely sorted.bam
rule map_stats:
    """Get the mapping stats for each aligment run"""
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        #uniq_bam="analysis/align/{sample}/{sample}_unique.sorted.bam"
    output:
        #temp("analysis/align/{sample}/{sample}_mapping.txt")
        "analysis/align/{sample}/{sample}_mapping.txt"
    threads: _align_threads
    message: "ALIGN: get mapping stats for each bam"
    log: _logfile
    #CAN/should samtools view be multi-threaded--
    #UPDATE: tricky on how to do this right w/ compounded commands
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
    log: _logfile
    run:
        files = " -f ".join(input)
        shell("cidc_wes/modules/scripts/align_getMapStats.py -f {files} > {output} 2>>{log}")

#REMOVING THIS
# rule sortUniqueBams:
#     """General sort rule--take a bam {filename}.bam and 
#     output {filename}.sorted.bam"""
#     input:
#         "analysis/align/{sample}/{sample}_unique.bam"
#     output:
#         #CANNOT temp this b/c it's used by qdnaseq!
#         "analysis/align/{sample}/{sample}_unique.sorted.bam",
#         #"analysis/align/{sample}/{sample}_unique.sorted.bam.bai"
#     message: "ALIGN: sort bam file"
#     log: _logfile
#     threads: _align_threads
#     shell:
#         "sambamba sort {input} -o {output} -t {threads} 2>>{log}"

#REPLACING picard with sentieon - next two rules
#NOTE: if we don't have a sentieon license, then should use sambamba markdup!!

rule scoreSample:
    "Calls sentieon driver  --fun score_info on the sample"
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    output:
        "analysis/align/{sample}/{sample}.sorted.score.txt"
    message: "ALIGN: score sample"
    log: _logfile
    threads: _align_threads
    params:
        index1=config['sentieon_path'],
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
    log: _logfile
    threads: _align_threads
    params:
        index1=config['sentieon_path'],
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {input.score} --metrics {output.met} {output.bamm}"""

# REMOVING this
# rule indexBam:
#     """Index bam file"""
#     input:
#         "analysis/align/{sample}/{prefix}_unique.{suffix}.bam"
#     output:
#         "analysis/align/{sample}/{prefix}_unique.{suffix}.bam.bai"
#     message: "ALIGN: indexing bam file {input}"
#     log: _logfile
#     threads: _align_threads
#     shell:
#         "sambamba index -t {threads} {input} {output}"

    
