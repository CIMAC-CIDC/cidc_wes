# module:Precision HLA typing from next-generation sequencing data by Optitype and Polysolver
_optitype_threads=16
# def polysolver_runsHelper(wildcards, iindex):
#     """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
#     returns the sample name of Normal (if iindex=0) else sample name of Tumor"""
#     tmp = []
#     r = config['runs'][wildcards.run]
#     #print(r)

#     #check that we have a valid pair
#     if len(r) >=2:
#         sample_name = r[iindex]
#         tmp.append("analysis/align/%s/%s.sorted.bam" % (sample_name, sample_name))
#     else:
#         #NOTE: I can't figure out a proper kill command so I'll do this
#         tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
#     #print(tmp)
#     return tmp

# def getNormal_sample(wildcards):
#     return polysolver_runsHelper(wildcards, 0)

# def getTumor_sample(wildcards):
#     return polysolver_runsHelper(wildcards, 1)

def optitype_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/optitype/%s/%s_result.tsv" % (sample,sample))
        ls.append("analysis/optitype/%s/%s_coverage_plot.pdf" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end1.fastq" % (sample,sample))
        ls.append("analysis/optitype/%s/%s.sorted.chr6.end2.fastq" % (sample,sample))
        # ls.append("analysis/optitype/%s/orig.winners.hla.txt" %(sample))
        # ls.append("analysis/optitype/%s/call_stats.$allele.out" %(sample))
        # ls.append("analysis/optitype/%s/$allele.all.somatic.indels.vcf" %(sample))
        # ls.append("analysis/optitype/%s/orig.indiv.mutect.unfiltered.nonsyn.annotated" %(sample))
        # ls.append("analysis/optitype/%s/orig.indiv.mutect.filtered.nonsyn.annotated" %(sample))
        # ls.append("analysis/optitype/%s/orig.indiv.mutect.syn.nonsyn.annotated" %(sample))
        # ls.append("analysis/optitype/%s/orig.indiv.mutect.ambiguous.annotated" %(sample))
        # ls.append("analysis/optitype/%s/orig.indiv.strelka_indels.filtered.annotated" %(sample))
        # ls.append("analysis/optitype/%s/orig.indiv.strelka_indels.ambiguous.annotated" %(sample))
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
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_bamtofastq.txt"
    shell:
        """bedtools bamtofastq -i {input.in_sortchr6bamfile} -fq {output.chr6fastqfile1} -fq2 {output.chr6fastqfile2}"""

rule optitype_hlatyping:
    """Precision HLA typing from next-generation sequencing data by OptiType
This will produce a time-stamped directory inside the specified outputn directory containing a CSV file with the predicted optimal (and if enumerated, sub-optimal)HLA genotype, and a pdf file containing a coverage plot of the predicted alleles for diagnostic purposes"""
    
    input:
        in_chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        in_chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    params:
        #PathtoOptiType = config['conda_path'],
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "./analysis/optitype/%s/" % (wildcards.sample),
        #outputname = lambda wildcards: wildcards.sample
        path="source activate /home/taing/miniconda3/envs/optitype/", #HARD-Coding the path and activateing conda env
        optitype_config="cidc_wes/static/optitype/config.ini",
    output:
        HLAgenotype = "analysis/optitype/{sample}/{sample}_result.tsv",
        Coverageplot = "analysis/optitype/{sample}/{sample}_coverage_plot.pdf"
    benchmark:
        "benchmarks/optitype/{sample}/{sample}.optitype_hlatyping.txt"
    shell:
        """{params.path}; OptiTypePipeline.py -i {input.in_chr6fastqfile1} {input.in_chr6fastqfile2} --dna -v -o {params.output_dir} -p {params.name} --config {params.optitype_config}"""
        
# rule polysolver:
#     input:
#         in_sortchr6bamfile = "analysis/HLATyping/{sample}/{sample}.sorted.chr6.bam"
#     params:
#         output_dir=lambda wildcards: "analysis/HLATyping/%s/" % (wildcards.sample)
#     singularity:
#         "docker://r/sachet/polysolver"
#     output:
#         HLAwinners = "analysis/HLATyping/{sample}/orig.winners.hla.txt"
#     shell:
#         """singularity exec polysolver.img /usr/local/libexec/polysolver/scripts/shell_call_hla_type -bam {input.in_sortchr6bamfile} Unknown 1 hg38 STDFQ 0 -outdir {params.output_dir}"""

# rule polysolver_mutation:
#     input:
#         tumor_bam=getTumor_sample,
#         normal_bam=getNormal_sample,
#         winners_file = "analysis/HLATyping/{sample}/orig.winners.hla.txt"
#     params:
#         output_dir=lambda wildcards: "analysis/HLATyping/%s/" % (wildcards.sample)
#     singularity:
#         "docker://r/sachet/polysolver"
#     output:
#         Mutect_output = "analysis/HLATyping/{sample}/call_stats.$allele.out",
#         Strelka_output = "analysis/HLATyping/{sample}/$allele.all.somatic.indels.vcf"
#     shell:
#         """singularity exec polysolver.img /usr/local/libexec/polysolver/scripts/shell_call_hla_mutations_from_type {input.normal_bam} {input.tumor_bam} {input.winners_file} -build hg38 -outDir {params.output_dir}"""

# rule polysolver_Annotation_mutation:
#     input:
#         input_filepath = "analysis/HLATyping/{sample}/"
#     singularity:
#         "docker://r/sachet/polysolver"
#     output:
#         Unfiltered.nonsyn.annotated = "analysis/HLATyping/{sample}/orig.indiv.mutect.unfiltered.nonsyn.annotated",
#         filetered.nonsyn.annotated = "analysis/HLATyping/{sample}/orig.indiv.mutect.filtered.nonsyn.annotated",
#         filtered.syn.annotated = "analysis/HLATyping/{sample}/orig.indiv.mutect.syn.nonsyn.annotated",
#         ambiguous.annotated = "analysis/HLATyping/{sample}/orig.indiv.mutect.ambiguous.annotated",
#         strelka_indels.filtered.annotated = "analysis/HLATyping/{sample}/orig.indiv.strelka_indels.filtered.annotated",
#         strelka_indels.ambiguous.annotated = "analysis/HLATyping/{sample}/orig.indiv.strelka_indels.ambiguous.annotated"
#     shell:
#         """singularity exec polysolver.img /usr/local/libexec/polysolver/scripts/shell_annotate_hla_mutations indiv {input.input_filepath}"""







