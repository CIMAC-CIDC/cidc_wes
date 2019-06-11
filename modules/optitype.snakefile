# module:Precision HLA typing from next-generation sequencing data by Optitype and Polysolver
_optitype_threads=16

def optitype_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/optitype/%s/%s_result.tsv" % (sample,sample))
        ls.append("analysis/optitype/%s/%s_coverage_plot.pdf" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end1.fastq" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end2.fastq" % (sample,sample))
    return ls
    
rule optitype_all:
    input:
        optitype_targets
        
rule optitype_extract_chr6:
    """Extract chr6 by sambamba"""
    input:
        in_sortbamfile = "analysis/align/{sample}/{sample}.sorted.bam"
    output:
        chr6sortbamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam"
    threads:_optitype_threads
    group: "optitype"
    conda: "../envs/optitype.yml"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_extract_chr6.txt"
    shell:
        """sambamba view -t {threads} -f bam -h {input.in_sortbamfile} chr6 > {output.chr6sortbamfile}"""
        
rule optitype_bamtofastq:
    """Convert the sorted.chr6.bam file to fastq by bedtools"""
    input:
        in_sortchr6bamfile = "analysis/optitype/{sample}/{sample}.sorted.chr6.bam"
    output:
        chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    group: "optitype"
    conda: "../envs/optitype.yml"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_bamtofastq.txt"
    shell:
        """bedtools bamtofastq -i {input.in_sortchr6bamfile} -fq {output.chr6fastqfile1} -fq2 {output.chr6fastqfile2}"""

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
        output_dir=lambda wildcards: "%sanalysis/optitype/%s/" % (config['remote_path'], wildcards.sample),
        #outputname = lambda wildcards: wildcards.sample
        path="source activate /home/taing/miniconda3/envs/optitype/", #HARD-Coding the path and activateing conda env
        optitype_config="cidc_wes/static/optitype/config.ini",
    output:
        HLAgenotype = "analysis/optitype/{sample}/{sample}_result.tsv",
        Coverageplot = "analysis/optitype/{sample}/{sample}_coverage_plot.pdf"
    group: "optitype"
    conda: "../envs/optitype.yml"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_hlatyping.txt"
    shell:
        """{params.path}; OptiTypePipeline.py -i {input.in_chr6fastqfile1} {input.in_chr6fastqfile2} --dna -v -o {params.output_dir} -p {params.name} --config {params.optitype_config}"""
