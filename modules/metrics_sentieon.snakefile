#MODULE: Data Metrics by Sentieon
#import os
_metrics_threads=8
##sentieon  config sentieon_path in the config file as below
#sentieon_path="/cluster/jxfu/proj/CIDC/Sentieon/release/sentieon-genomics-201808.01/bin/"
#export SENTIEON_LICENSE=172.24.216.24:8990

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
	
    #TODO- fill this in
    return ls

rule metrics_all:
    input:
        metrics_targets

rule Metrics_sentieon:
    """Get the metrics calculations from  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.bam",
         bai="analysis/align/{sample}/{sample}.sorted.bam.bai"
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
    shell:
        """{params.index1}/sentieon plot metrics -o {output.overallmetricspdf}  gc={input.gcmetrics} qd={input.qd} mq={input.mq}  isize={input.insertsize}"""
