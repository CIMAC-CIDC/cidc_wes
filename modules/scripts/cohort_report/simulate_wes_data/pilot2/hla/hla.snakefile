def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        tmr = config[run]['tumor']
        nrm = config[run]['normal']
        ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (tmr,tmr))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (nrm,nrm))
        
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end1.fastq" % (tmr,tmr))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end2.fastq" % (tmr,tmr))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end1.fastq" % (nrm,nrm))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end2.fastq" % (nrm,nrm))

        ls.append("analysis/optitype/%s/%s_result.tsv" % (tmr,tmr))
        ls.append("analysis/optitype/%s/%s_result.tsv" % (nrm,nrm))

        ls.append("analysis/xhla/%s/report-%s-hla.json" % (tmr,tmr))
        ls.append("analysis/xhla/%s/report-%s-hla.json" % (nrm,nrm))
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule optitype_extract_chr6:
    """Extract chr6 by sambamba"""
    input:
        in_sortbamfile = "raw_data/{sample}.sorted.dedup.bam"
    output:
        chr6sortbamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam"
    threads: 8
    group: "optitype"
    #conda: "../envs/optitype.yml"
    #benchmark:
    #    "benchmarks/optitype/{sample}/{sample}.optitype_extract_chr6.txt"
    shell:
        """sambamba view -t {threads} -f bam -h {input.in_sortbamfile} chr6 > {output.chr6sortbamfile}"""

rule optitype_index_chr6bam:
    """index chr6bam"""
    input:
        "analysis/optitype/{sample}/{sample}.sorted.chr6.bam"
    output:
        "analysis/optitype/{sample}/{sample}.sorted.chr6.bam.bai"
    threads: 8
    group: "optitype"
    #conda: "../envs/optitype.yml"
    #benchmark:
    #    "benchmarks/optitype/{sample}/{sample}.optitype_index_chr6bam.txt"
    shell:
        """sambamba index -t {threads} {input}"""

rule optitype_bamtofastq:
    """Convert the sorted.chr6.bam file to fastq by samtools"""
    input:
        in_sortchr6bamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam",
        in_index = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam.bai"
    output:
        chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    group: "optitype"
    #conda: "../envs/optitype.yml"
    log: "analysis/logs/optitype/{sample}/{sample}.optitype_bamtofastq.log"
    #benchmark:
    #    "benchmarks/optitype/{sample}/{sample}.optitype_bamtofastq.txt"
    shell:
        """samtools fastq -@ 2 -1 {output.chr6fastqfile1} -2 {output.chr6fastqfile2} {input.in_sortchr6bamfile} 2> {log}"""

rule optitype_hlatyping:
    """Precision HLA typing from next-generation sequencing data by
    OptiType This will produce a time-stamped directory inside the
    specified outputn directory containing a CSV file with the predicted
    optimal (and if enumerated, sub-optimal)HLA genotype, and a pdf file
    containing a coverage plot of the predicted alleles for diagnostic
    purposes"""
    
    input:
        in_chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        in_chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    params:
        #PathtoOptiType = config['conda_path'],
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "analysis/optitype/%s/" % wildcards.sample,
        #outputname = lambda wildcards: wildcards.sample
        #path="source activate %s" % config['optitype_root'],
        path="source activate ~/miniconda3/envs/optitype",
        optitype_config="/mnt/ssd/mda-r1-pt1_report/cidc_wes/static/optitype/config.ini",
    output:
        HLAgenotype = "analysis/optitype/{sample}/{sample}_result.tsv",
        Coverageplot = "analysis/optitype/{sample}/{sample}_coverage_plot.pdf"
    group: "optitype"
    #conda: "../envs/optitype.yml"
    #benchmark:
    #    "benchmarks/optitype/{sample}/{sample}.optitype_hlatyping.txt"
    shell:
        """{params.path}; OptiTypePipeline.py -i {input.in_chr6fastqfile1} {input.in_chr6fastqfile2} --dna -v -o {params.output_dir} -p {params.name} --config {params.optitype_config}"""

rule xhla:
    """calculate hlatyping by xhla"""
    input:
        #in_sortbamfile = "analysis/align/{sample}/{sample}.sorted.dedup.bam"
        #LEN- changing this to use chrom6 as well
        in_sortbamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam",
        in_index = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam.bai",
    output:
        chr6sortbamfile = "analysis/xhla/{sample}/report-{sample}-hla.json"
    threads: 8
    group: "xhla"
    params:
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "analysis/xhla/%s/" % wildcards.sample,
        #path="source activate %s" % config['xhla_root'],
        path="source activate ~/miniconda3/envs/xHLA",
    #singularity: "docker://humanlongevity/hla"
    #conda: "../envs/xhla_env.yml"
    #benchmark:
    #    "benchmarks/xhla/{sample}/{sample}.xhla.txt"
    shell:
        """{params.path}; run.py --sample_id {params.name} --input_bam_path {input.in_sortbamfile} --output_path  {params.output_dir}"""
