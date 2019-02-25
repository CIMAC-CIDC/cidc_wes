#MODULE: Coverage Metrics by Sentieon
#import os
_coveragemetrics_threads=16


def coveragemetrics_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/metrics/%s/%s_coverage_metrics.txt" % (sample,sample))
    return ls

rule coveragemetrics_all:
    input:
        coveragemetrics_targets

rule CoverageMetrics_sentieon:
    """Get the metrics calculations from  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.bam",
         bai="analysis/align/{sample}/{sample}.sorted.bam.bai",
    output:
         coveragemetrics="analysis/metrics/{sample}/{sample}_coverage_metrics.txt",
    message:
        "COVERAGE: coverage calculations for bam file"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        cov_thresh=50, #LT: put this in config.yaml
        index2=config['CDS_bed_input'],
    threads: _coveragemetrics_threads
    benchmark:
        "benchmarks/coverage/{sample}/{sample}.CoverageMetrics.txt"
    shell:
        #change cov_thresh as an user input config
        """{params.index1}/sentieon driver -r {params.index}  -t  {threads} --interval {params.index2} -i {input.bam} --algo CoverageMetrics --cov_thresh {params.cov_thresh} {output.coveragemetrics}"""
