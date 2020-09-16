#MODULE: Data Metrics by Sentieon
#import os
_metrics_threads=32

def metrics_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/metrics/%s/%s_mq_metrics.txt" % (sample,sample))
	ls.append("analysis/metrics/%s/%s_qd_metrics.txt" % (sample,sample))
	ls.append("analysis/metrics/%s/%s_gc_summary.txt" % (sample,sample))
	ls.append("analysis/metrics/%s/%s_gc_metrics.txt" % (sample,sample))
	ls.append("analysis/metrics/%s/%s_aln_metrics.txt" % (sample,sample))
	ls.append("analysis/metrics/%s/%s_is_metrics.txt" % (sample,sample))
	ls.append("analysis/metrics/%s/%s_metrics.pdf" % (sample,sample))

        #JSON files
	ls.append("analysis/report/json/gc_content/%s.gc.json" % sample)
        ls.append("analysis/report/json/insert_size/%s.insert_size.json" % sample)
    #coverage metrics summaries
    ls.append("analysis/metrics/all_sample_summaries.txt")
    return ls

rule metrics_all:
    input:
        metrics_targets

rule Metrics_sentieon:
    """Get the metrics calculations from  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.dedup.bam",
         bai="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai"
    output:
         mq="analysis/metrics/{sample}/{sample}_mq_metrics.txt",
         qd="analysis/metrics/{sample}/{sample}_qd_metrics.txt",
         gcsummary="analysis/metrics/{sample}/{sample}_gc_summary.txt",
         gcmetrics="analysis/metrics/{sample}/{sample}_gc_metrics.txt",
         aln="analysis/metrics/{sample}/{sample}_aln_metrics.txt",
         insertsize="analysis/metrics/{sample}/{sample}_is_metrics.txt",
    message:
        "METRICS: quality control for  mapped reads"
    #log: _logfile
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
    threads: _metrics_threads
    group: "metrics"
    benchmark:
        "benchmarks/metrics/{sample}/{sample}.Metrics_sentieon.txt"
    shell:
        """{params.index1}/sentieon driver -r {params.index}  -t  {threads} -i {input.bam} --algo MeanQualityByCycle {output.mq} --algo QualDistribution {output.qd}  --algo GCBias --summary {output.gcsummary}  {output.gcmetrics} --algo AlignmentStat --adapter_seq '' {output.aln} --algo InsertSizeMetricAlgo {output.insertsize}"""


rule Metrics_sentieon_plots:
    """Get the metrics plots"""
    input:
        mq="analysis/metrics/{sample}/{sample}_mq_metrics.txt",
        qd="analysis/metrics/{sample}/{sample}_qd_metrics.txt",
        gcsummary="analysis/metrics/{sample}/{sample}_gc_summary.txt",
        gcmetrics="analysis/metrics/{sample}/{sample}_gc_metrics.txt",
        aln="analysis/metrics/{sample}/{sample}_aln_metrics.txt",
        insertsize="analysis/metrics/{sample}/{sample}_is_metrics.txt",
    output:
        overallmetricspdf="analysis/metrics/{sample}/{sample}_metrics.pdf"
    message:
        "METRICS: plot data metrics for mapped reads"
    #log: _logfile
    params:
        index1=config['sentieon_path'],
    threads: _metrics_threads
    group: "metrics"
    benchmark:
        "benchmarks/metrics/{sample}/{sample}.Metrics_sentieon_plots.txt"
    shell:
        """{params.index1}/sentieon plot metrics -o {output.overallmetricspdf}  gc={input.gcmetrics} qd={input.qd} mq={input.mq}  isize={input.insertsize}"""

rule metrics_collect_target_summaries:
    """Collect all of the sample summaries and put them in one file
    input: {sample}_target_metrics.sample_summary.txt from coverage.snakefile
    """
    input:
        coverage=expand("analysis/metrics/{sample}/{sample}_target_metrics.sample_summary.txt", sample=sorted(config['samples'])),
        align=expand("analysis/align/{sample}/{sample}_mapping.txt", sample=sorted(config['samples'])),
    output:
        "analysis/metrics/all_sample_summaries.txt"
    params:
        cov_files = lambda wildcards, input: " -f ".join(input.coverage),
        align_files = lambda wildcards, input: " -a ".join(input.align)
    group: "metrics"
    benchmark:
        "benchmarks/metrics/metrics_collect_target_summaries.txt"
    shell:
        "cidc_wes/modules/scripts/metrics_collect_target_summaries.py -f {params.cov_files} -a {params.align_files} > {output}"

rule metrics_json_gc_content:
    """jsonify the GC content contained in {sample}/{sample}_gc_metrics.txt
    """
    input:
        "analysis/metrics/{sample}/{sample}_gc_metrics.txt",
    output:
        "analysis/report/json/gc_content/{sample}.gc.json"
    group: "metrics"
    benchmark:
        "benchmarks/metrics/{sample}.metrics_json_gc_content.txt"
    shell:
        "cidc_wes/modules/scripts/json_gc_content.py -f {input} -o {output}"

rule metrics_json_insert_size:
    """jsonify the insert size contained in {sample}/{sample}_is_metrics.txt
    """
    input:
        "analysis/metrics/{sample}/{sample}_is_metrics.txt",
    output:
        "analysis/report/json/insert_size/{sample}.insert_size.json"
    group: "metrics"
    benchmark:
        "benchmarks/metrics/{sample}.metrics_json_insert_size.txt"
    shell:
        "cidc_wes/modules/scripts/json_insert_size.py -f {input} -o {output}"
