#MODULE: Coverage Metrics by Sentieon
#import os
_coveragemetrics_threads=8


def coveragemetrics_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/metrics/%s/%s_coverage_metrics.txt" % (sample,sample))

    #TODO- fill this in
    return ls

rule metrics_all:
    input:
        coveragemetrics_targets

rule CoverageMetrics_sentieon:
    """Get the metrics calculations from  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.bam"
    output:
         coveragemetrics="analysis/metrics/{sample}/{sample}_coverage_metrics.txt",
    message:
        "COVERAGE: coverage calculations for bam file"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
    threads: _coveragemetrics_threads
    shell:
        """{params.index1}/sentieon driver -r {params.index}  -t  {threads} -i {input.bam} --algo CoverageMetrics --cov_thresh 50 {output.coveragemetrics}"""
