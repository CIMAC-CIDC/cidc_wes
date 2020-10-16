#_base_dir="/mnt/ssd/mda-r1-pt1_report/"
_base_dir="/mnt/ssd/mda-r1-pt1_report/simulate_data/pilot2/coverage_full/"

def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        tmr = config[run]['tumor']
        nrm = config[run]['normal']
        ls.append("analysis/metrics/%s/%s_coverage_metrics.sample_summary.txt" % (tmr,tmr))
        ls.append("analysis/metrics/%s/%s_coverage_metrics.sample_summary.txt" % (nrm,nrm))

        ls.append("json/%s.coverage.json" % run)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())


configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule CoverageMetrics_sentieon:
    """Get the metrics calculations from  mapped reads"""
    input:
         bam="raw_data/{sample}.sorted.dedup.bam",
         bai="raw_data/{sample}.sorted.dedup.bam.bai",
    output:
         coveragemetrics="analysis/metrics/{sample}/{sample}_coverage_metrics.txt",
         coveragemetrics_summary="analysis/metrics/{sample}/{sample}_coverage_metrics.txt.sample_summary",
    message:
        "COVERAGE: coverage calculations for bam file"
    params:
        #index=config['genome_fasta'],
        #index= os.path.join(_base_dir, "ref_files/hg38/bwa_gdc/GRCh38.d1.vd1.CIDC.fa"),
        index= os.path.join(_base_dir, "ref_files/hg38/bwa_gdc/GRCh38.d1.vd1.fa"),
        #index1=config['sentieon_path'],
        index1="/home/taing/sentieon/sentieon-genomics-201808.05/bin/",
        cov_thresh=50, #LT: put this in config.yaml
        #index2=config['CDS_Bed_input'],
        index2=os.path.join(_base_dir, "./ref_files/hg38/gencode27.canonical.bed"),
    threads: 4 #16
    benchmark:
        "benchmarks/coverage/{sample}/{sample}.CoverageMetrics.txt"
    shell:
        #change cov_thresh as an user input config
        """{params.index1}/sentieon driver -r {params.index}  -t  {threads} --interval {params.index2} -i {input.bam} --algo CoverageMetrics --cov_thresh {params.cov_thresh} {output.coveragemetrics}"""

rule addExtension:
    """Simple rule to add .txt to our file"""
    input:
        "analysis/metrics/{sample}/{sample}_coverage_metrics.txt.sample_summary",
    output:
        "analysis/metrics/{sample}/{sample}_coverage_metrics.sample_summary.txt"
    shell:
        "mv {input} {output}"

def json_getFiles(wildcards):
    """Will return a list of the run's tumor and normal files"""
    tmp = {}
    run = config[wildcards.run]
    tmr = run['tumor']
    nrm = run['normal']
    
    tmp['tumor'] = "analysis/metrics/%s/%s_coverage_metrics.sample_summary.txt" % (tmr, tmr)
    tmp['normal'] = "analysis/metrics/%s/%s_coverage_metrics.sample_summary.txt" % (nrm, nrm)

    #print(tmp)
    return tmp

rule coverage_content_json:
    input:
        unpack(json_getFiles)
    output:
        "json/{run}.coverage.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./coverage_json_writer.py -r {params.run} -t {input.tumor} -n {input.normal} -o {output}"
