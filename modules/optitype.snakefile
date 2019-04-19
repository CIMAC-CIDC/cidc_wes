# module:Precision HLA typing from next-generation sequencing data by Optitype
sample = "mocha2-Run2-pt7-Tumor"
def optitype_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/HLATyping/%s/%s_result.tsv" % (sample,sample))
        ls.append("analysis/HLATyping/%s/%s_coverage_plot.pdf" % (sample,sample))
        ls.append("analysis/HLATyping/%s/%s.sorted.chr6.bam" % (sample,sample))
        ls.append("analysis/HLATyping/%s/%s.sorted.chr6.end1.fastq" % (sample,sample))
        ls.append("analysis/HLATyping/%s/%s.sorted.chr6.end2.fastq" % (sample,sample))
    return ls
    
rule optitype_all:
    input:
        optitype_targets
        
rule optitype_extract_chr6:
    """Extract chr6 by samtool"""
    input:
        in_sortbamfile = "analysis/align/{sample}/{sample}.sorted.bam"
    output:
        chr6sortbamfile = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.bam"
    shell:
        """samtools view -b -h {input.in_sortbamfile} chr6 > {output.chr6sortbamfile}"""
        
rule optitype_bamtofastq:
    """Convert the sorted.chr6.bam file to fastq by bedtools"""
    input:
        in_sortchr6bamfile = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.bam"
    output:
        chr6fastqfile1 = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.end1.fastq",
        chr6fastqfile2 = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.end2.fastq"
    shell:
        """bedtools bamtofastq -i {input.in_sortchr6bamfile} -fq {output.chr6fastqfile1} -fq2 {output.chr6fastqfile2}"""

rule optitype_hlatyping:
    """Precision HLA typing from next-generation sequencing data by OptiType
This will produce a time-stamped directory inside the specified outputn directory containing a CSV file with the predicted optimal (and if enumerated, sub-optimal)HLA genotype, and a pdf file containing a coverage plot of the predicted alleles for diagnostic purposes"""
    
    input:
        in_chr6fastqfile1 = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.end1.fastq",
        in_chr6fastqfile2 = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.end2.fastq"
    params:
        #PathtoOptiType = config['optitype_path'],
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "./analysis/HLATyping/%s/" % (wildcards.sample)
        #outputname = lambda wildcards: wildcards.sample
    output:
        HLAgenotype = "analysis/HLATyping/{sample}/{sample}_result.tsv",
        Coverageplot = "analysis/HLATyping/{sample}/{sample}_coverage_plot.pdf"
    shell:
        """python /homes/jjyu/miniconda3/bin/OptiTypePipeline.py -i {input.in_chr6fastqfile1} {input.in_chr6fastqfile2} --dna -v -o {params.output_dir} -p {params.name} """
