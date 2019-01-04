#MODULE: Align fastq files to genome - common rules
#import os
_align_threads=8

def align_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/align/%s/%s.sorted.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s.sorted.bam.bai" % (sample,sample))
        ls.append("analysis/align/%s/%s_unique.sorted.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_unique.sorted.bam.bai"%(sample,sample))
        ls.append("analysis/align/%s/%s_unique.sorted.dedup.bam" % (sample,sample))
        ls.append("analysis/align/%s/%s_unique.sorted.dedup.bam.bai" % (sample,sample))
        #ls.append("analysis/align/%s/%s.unmapped.fq.gz" % (sample,sample))
        ls.append("analysis/align/%s/%s_readsPerChrom.txt" % (sample,sample))
    ls.append("analysis/align/mapping.csv")
    return ls


def getBam(wildcards):
    """This input fn will check to see if the user specified a .fastq or a .bam
    for the sample.  IF the former (.fastq), will simply return the canonical
    path, otherwise (.bam) will return the user-specified (bam) path"""
    #CHECK first entry's file suffix
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    ret = "analysis/align/%s/%s.bam" % (s,s)
    if first_file.endswith('.bam'):
        #CLEANER to check for .bam vs (.fastq, fq, fq.gz, etc)
        ret = first_file
    return [ret]

rule align_all:
    input:
        align_targets

rule uniquely_mapped_reads:
    """Get the uniquely mapped reads"""
    input:
        #"analysis/align/{sample}/{sample}.bam"
        getBam
    output:
        temp("analysis/align/{sample}/{sample}_unique.bam")
    message: "ALIGN: Filtering for uniquely mapped reads"
    log: _logfile
    threads: _align_threads
    shell:
        #NOTE: this is the generally accepted way of doing this as multiply 
        #mapped reads have a Quality score of 0
        #NOTE: -@ = --threads
        "samtools view -bq 1 -@ {threads} {input} > {output}"

rule map_stats:
    """Get the mapping stats for each aligment run"""
    input:
        #bam="analysis/align/{sample}/{sample}.bam",
        bam=getBam,
        uniq_bam="analysis/align/{sample}/{sample}_unique.bam"
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
        "samtools flagstat {input.bam} > {output} 2>>{log}"
        " && samtools view -c {input.uniq_bam} >> {output} 2>> {log}"

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

rule sortBams:
    """General sort rule--take a bam {filename}.bam and 
    output {filename}.sorted.bam"""
    input:
        #"analysis/align/{sample}/{filename}.bam"
        getBam
    output:
        "analysis/align/{sample}/{sample}.sorted.bam",
        #"analysis/align/{sample}/{sample}.sorted.bam.bai"
    message: "ALIGN: sort bam file"
    log: _logfile
    threads: _align_threads
    shell:
        "sambamba sort {input} -o {output} -t {threads} 2>>{log}"

rule sortUniqueBams:
    """General sort rule--take a bam {filename}.bam and 
    output {filename}.sorted.bam"""
    input:
        "analysis/align/{sample}/{sample}_unique.bam"
    output:
        #CANNOT temp this b/c it's used by qdnaseq!
        "analysis/align/{sample}/{sample}_unique.sorted.bam",
        #"analysis/align/{sample}/{sample}_unique.sorted.bam.bai"
    message: "ALIGN: sort bam file"
    log: _logfile
    threads: _align_threads
    shell:
        "sambamba sort {input} -o {output} -t {threads} 2>>{log}"

rule dedupSortedUniqueBams:
    """Dedup sorted unique bams using PICARD
    output {sample}_unique.sorted.dedup.bam"""
    input:
        "analysis/align/{sample}/{sample}_unique.sorted.bam"
    output:
        "analysis/align/{sample}/{sample}_unique.sorted.dedup.bam"
    message: "ALIGN: dedup sorted unique bam file"
    log: _logfile
    threads: _align_threads
    shell:
        "picard MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE={log} 2>> {log}"

rule indexBam:
    """Index bam file"""
    input:
        "analysis/align/{sample}/{prefix}.bam"
    output:
        "analysis/align/{sample}/{prefix}.bam.bai"
    message: "ALIGN: indexing bam file {input}"
    log: _logfile
    shell:
        "sambamba index {input} {output}"

#SKIP- DEPRECATED
# rule extractUnmapped:
#     """Extract the unmapped reads and save as {sample}.unmapped.bam"""
#     input:
#         #"analysis/align/{sample}/{sample}.bam"
#         getBam
#     output:
#         temp("analysis/align/{sample}/{sample}.unmapped.bam")
#     message: "ALIGN: extract unmapped reads"
#     log: _logfile
#     threads: _align_threads
#     shell:
#         #THIS extracts all unmapped reads
#         #"samtools view -b -f 4 --threads {threads} {input} >{output} 2>>{log}"
#         #THIS extracts all READ (pairs) where at least one in unmapped
#         #ref: https://www.biostars.org/p/56246/ search "rgiannico"
#         #NOTE: -@ = --threads
#         "samtools view -b -F 2 -@ {threads} {input} > {output} 2>>{log}"

# rule bamToFastq:
#     """Convert the unmapped.bam to fastq"""
#     input:
#         "analysis/align/{sample}/{sample}.unmapped.bam"
#     output:
#         "analysis/align/{sample}/{sample}.unmapped.fq"
#     params:
#         #handle PE alignments!
#         mate2 = lambda wildcards: "-fq2 analysis/align/{sample}/{sample}.unmapped.fq2" if len(config["samples"][wildcards.sample]) == 2 else ""
#     message: "ALIGN: convert unmapped bam to fastq"
#     log: _logfile
#     shell:
#         "bamToFastq -i {input} -fq {output} {params.mate2}"

# rule gzipUnmappedFq:
#     """gzip unmapped fq(s)"""
#     input:
#         "analysis/align/{sample}/{sample}.unmapped.fq"
#     output:
#         "analysis/align/{sample}/{sample}.unmapped.fq.gz"
#     params:
#         #handle PE alignments!
#         mate2 = lambda wildcards: "analysis/align/{sample}/{sample}.unmapped.fq2" if len(config["samples"][wildcards.sample]) == 2 else ""
#     message: "ALIGN: gzip unmapped fq files"
#     log: _logfile
#     shell:
#         "gzip {input} {params} 2>>{log}"

rule readsPerChromStat:
    """For each sample, generates a _readsPerChrom.txt file, which is:
    chr1   #readsOnChr1
    ...
    chrX   #readsOnChrX
    """
    input:
        bam = "analysis/align/{sample}/{sample}.sorted.bam",
        #NOTE: even though we don't use the bai, we need to ensure bam sorted
        bai = "analysis/align/{sample}/{sample}.sorted.bam.bai"
    params:
        awk_call = """awk '{print $1 \"\\t\" $3; s+=$3} END {print \"total reads = \" s}'"""
    output:
        "analysis/align/{sample}/{sample}_readsPerChrom.txt"
    message: "ALIGN: collecting the number of reads per chrom"
    log: _logfile
    shell:
        "cidc_wes/modules/scripts/align_readsPerChrom.sh -a {input.bam} > {output} 2>> {log}"

    
