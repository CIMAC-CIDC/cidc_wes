#MODULE: Coverage Metrics by Sentieon
#import os
_coveragemetrics_threads=16


def coveragemetrics_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
    	ls.append("analysis/metrics/%s/%s_coverage_metrics.txt" % (sample,sample))
        ls.append("analysis/metrics/%s/%s_target_metrics.txt" % (sample,sample))
        for center in center_targets:
            ls.append("analysis/metrics/%s/%s.%s.mosdepth.region.dist.txt" % (sample,sample,center))
            ls.append("analysis/metrics/%s/%s.%s.mosdepth.region.summary.txt" % (sample,sample,center))
    return ls

rule coverage_all:
    input:
        coveragemetrics_targets

rule CoverageMetrics_sentieon:
    """Get the metrics calculations from  mapped reads"""
    input:
         bam="analysis/align/{sample}/{sample}.sorted.dedup.bam",
         bai="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
    output:
         coveragemetrics="analysis/metrics/{sample}/{sample}_coverage_metrics.txt",
    message:
        "COVERAGE: coverage calculations for bam file"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        cov_thresh=50, #LT: put this in config.yaml
        index2=config['CDS_Bed_input'],
    threads: _coveragemetrics_threads
    group: "coverage"
    benchmark:
        "benchmarks/coverage/{sample}/{sample}.CoverageMetrics.txt"
    shell:
        #change cov_thresh as an user input config
        """{params.index1}/sentieon driver -r {params.index}  -t  {threads} --interval {params.index2} -i {input.bam} --algo CoverageMetrics --cov_thresh {params.cov_thresh} {output.coveragemetrics}"""


rule targets_sentieon:
    """Get the coverage metrics from target beds"""
    input:
        bam="analysis/align/{sample}/{sample}.sorted.dedup.bam",
        bai="analysis/align/{sample}/{sample}.sorted.dedup.bam.bai",
    output:
        targetmetrics="analysis/metrics/{sample}/{sample}_target_metrics.txt",
        summary="analysis/metrics/{sample}/{sample}_target_metrics.txt.sample_summary",
    message:
        "Coverage calculation from target bed files"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        cov_thresh=50,
        index2=config['target_Bed_input'],
    threads: _coveragemetrics_threads
    group: "coverage"
    benchmark:
        "benchmarks/targetcoverage/{sample}/{sample}.targetMetrics.txt"
    shell:
        """{params.index1}/sentieon driver -r {params.index}  -t  {threads} --interval {params.index2} -i {input.bam} --algo CoverageMetrics --cov_thresh {params.cov_thresh} {output.targetmetrics}"""


rule coverage_mosdepth:
    """SAMPLE coverage at target sites"""
    input:
        bam="analysis/align/{sample}/{sample}_recalibrated.bam",
        bai="analysis/align/{sample}/{sample}_recalibrated.bam.bai"
    params:
        #FOUND center_targets found in wes.snakefile
        target= lambda wildcards: center_targets[wildcards.center],
        prefix=lambda wildcards: "%sanalysis/metrics/%s/%s.%s" % (config['remote_path'], wildcards.sample, wildcards.sample, wildcards.center),
    output:
        "analysis/metrics/{sample}/{sample}.{center}.mosdepth.region.dist.txt",
    group: "coverage"
    conda: "../envs/coverage.yml"
    threads: _coveragemetrics_threads
    benchmark:
        "benchmarks/metrics/{sample}/{sample}.{center}.mosdepth.txt"
    shell:
        "mosdepth -t {threads} -b {params.target} {params.prefix} {input.bam}"

rule coverage_summarize_sample_depth:
    """Summarize the results of the coverage_mosdepth output, 
    e.g. % 20x, 50x, etc"""
    input:
        "analysis/metrics/{sample}/{sample}.{center}.mosdepth.region.dist.txt"
    output:
        "analysis/metrics/{sample}/{sample}.{center}.mosdepth.region.summary.txt"
    params:
        sample_name=lambda wildcards:"%s.%s" % (wildcards.sample, wildcards.center)
    group: "coverage"
    conda: "../envs/coverage.yml"
    benchmark:
        "benchmarks/metrics/{sample}/{sample}.{center}.coverage_summarize_sample_depth.txt"
    shell:
        "Rscript cidc_wes/modules/scripts/mosdepth_summary.R {input} {params.sample_name} {output}"
    
