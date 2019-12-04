#MODULE: Align fastq files to genome - common rules
#import os
#_logfile="analysis/logs/align.log"
_align_threads=32
_bwa_threads=32

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

def align_getBam(wildcards):
    bam = config["samples"][wildcards.sample][0] #returns only the first elm
    return bam

rule align_all:
    input:
        align_targets

rule map_all:
    input:
        align_mapping_targets

#BASED on xindong's aggregate_align_input in CIDC_Chips
def aggregate_align_input(wildcards):
    # handle .bam files separately from .fastq files
    #check only the first file
    sample_first_file = config["samples"][wildcards.sample][0]
    if sample_first_file.endswith(".bam"):
        return ["analysis/align/{sample}/{sample}.sorted.fromBam.bam",
                "analysis/align/{sample}/{sample}.sorted.fromBam.bam.bai"]
    else:
        return ["analysis/align/{sample}/{sample}.sorted.fromFastq.bam",
                "analysis/align/{sample}/{sample}.sorted.fromFastq.bam.bai"]

rule aggregate_input:
    input:
        aggregate_align_input
    params:
        bam = lambda wildcards,input: input[0],
        bai = lambda wildcards,input: input[1],
    output:
       bam="analysis/align/{sample}/{sample}.sorted.bam",
       bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    shell:
        "mv {params.bam} {output.bam} && mv {params.bai} {output.bai}"

rule align_from_bam:
    input:
        align_getBam
    output:
        bam="analysis/align/{sample}/{sample}.sorted.fromBam.bam",
        bai="analysis/align/{sample}/{sample}.sorted.fromBam.bam.bai"
    threads: 32
    params:
        sentieon_path=config['sentieon_path'],
        bwa_index=config['bwa_index'],
        #DON'T write a mini-program withi a program-
        #awk cmd to add sample names to RGs!!
        awk_cmd=lambda wildcards: "awk -v OFS=\'\\t\' \'{ split($2,a,\":\"); read_id=a[2]; $2=\"ID:%s.\" read_id; print $0}\'" % wildcards.sample,
        #NEVER do it twice!- gawk cmd to inject sample name into each read!!!
        gawk_cmd=lambda wildcards: "gawk -v OFS=\'\\t\' \'{rg=match($0,/RG:Z:(\S+)/,a); read_id=a[1]; if (rg) {sub(/RG:Z:\S+/, \"RG:Z:%s.\" read_id, $0); print $0} else { print $0 }}\'" % wildcards.sample,
    benchmark: "benchmarks/align/{sample}/{sample}.align_from_bam.txt"
    shell:
        """samtools view -H {input} | grep \"^@RG\" | {params.awk_cmd} > {wildcards.sample}.header && samtools collate --output-fmt SAM -@ {threads} -Of {input} | {params.gawk_cmd} | samtools view -@ {threads} -b - | samtools fastq -@ {threads} -t -s /dev/null -0 /dev/null - | ({params.sentieon_path}/sentieon bwa mem -t {threads} -M -K 10000000 -p -C -H {wildcards.sample}.header {params.bwa_index} - || echo -n 'error' ) | {params.sentieon_path}/sentieon util sort -t {threads} -o {output.bam} --sam2bam -; rm {wildcards.sample}.header"""

rule align_from_fastq:
    input:
        align_getFastq
    output:
        bam="analysis/align/{sample}/{sample}.sorted.fromFastq.bam",
        bai="analysis/align/{sample}/{sample}.sorted.fromFastq.bam.bai"
    params:
        sentieon_path=config['sentieon_path'],
        bwa_index=config['bwa_index'],
        read_group= lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample),        
        input_bases="10000000",
        #need to adjust threads for the other process
        tthreads=lambda wildcards, input, output, threads, resources: threads-1
    threads: _bwa_threads
    message: "ALIGN: Running sentieon BWA mem for alignment"
    log: "analysis/logs/align/{sample}/align.sentieon_bwa.{sample}.log"
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.align_from_fastq.txt"
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
    log: "analysis/logs/align/{sample}/align.map_stats.{sample}.log"
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
    log: "analysis/logs/align/align.collect_map_stats.log"
    params:
        files = lambda wildcards, input: " -f ".join(input)
    group: "align"
    benchmark:
        "benchmarks/align/align_collect_map_stats.txt"
    shell:
        "cidc_wes/modules/scripts/align_getMapStats.py -f {params.files} > {output} 2>>{log}"
    #run:
    #    files = " -f ".join(input)
    #    shell("cidc_wes/modules/scripts/align_getMapStats.py -f {files} > {output} 2>>{log}")

rule scoreSample:
    "Calls sentieon driver  --fun score_info on the sample"
    input:
        bam="analysis/align/{sample}/{sample}.sorted.bam",
        bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    output:
        score="analysis/align/{sample}/{sample}.sorted.score.txt",
        idx="analysis/align/{sample}/{sample}.sorted.score.txt.idx",
    message: "ALIGN: score sample"
    log: "analysis/logs/align/{sample}/align.scoreSample.{sample}.log"
    threads: _align_threads
    params:
        index1=config['sentieon_path'],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.scoreSample.txt"
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.score}"""

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
    log: "analysis/logs/align/{sample}/align.dedupSortedUniqueBam.{sample}.log"
    threads: _align_threads
    params:
        index1=config['sentieon_path'],
    group: "align"
    benchmark:
        "benchmarks/align/{sample}/{sample}.dedupSortedUniqueBam.txt"
    shell:
        """{params.index1}/sentieon driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {input.score} --metrics {output.met} {output.bamm}"""
