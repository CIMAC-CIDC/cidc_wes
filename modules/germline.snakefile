# module: Generate SNP92 to identification of germline variant.
# paper_name:SNPs for a universal individual identification panel

def germline_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        #Consolidate these with an inner-for-loop?
        ls.append("analysis/germline/%s/%s_variant.vcf" % (sample,sample))
        ls.append("analysis/germline/%s/%s_SNP92.recode.vcf" % (sample,sample))
    return ls

rule germline_all:
    input:
        germline_targets

rule germline_bamtovcf:
    input:
        input_sortbamfile = "analysis/align/{sample}/{sample}.sorted.bam" 
    output:
        output_vcf="analysis/germline/{sample}/{sample}_variant.vcf"
    params:
        index=config['genome_fasta'],
        positions_bamtovcf=config['positions_bamtovcf']
    group: "germline"
    conda: 
    	"../envs/germline.yml"
    benchmark:
        "benchmarks/germline/{sample}/{sample}.bamtovcf.txt"
    shell:
        """samtools mpileup -g -Q 0 -f {params.index} {input.input_sortbamfile} --positions {params.positions_bamtovcf} | bcftools view > {output.output_vcf}"""

rule germline_snp92:
    input:
        input_vcf="analysis/germline/{sample}/{sample}_variant.vcf"
    output:
        output_SNP92="analysis/germline/{sample}/{sample}_SNP92.recode.vcf"
    params:
        positons_SNP92=config['positions_SNP92'],
        outname=lambda wildcards: "%sanalysis/germline/%s/%s_SNP92" % (config['remote_path'], wildcards.sample, wildcards.sample),
    group: "germline"
    conda:
        "../envs/germline.yml"
    benchmark:
        "benchmarks/germline/{sample}/{sample}.snp92.txt"
    shell:
        """vcftools --vcf {input.input_vcf} --positions {params.positons_SNP92} --recode --out {params.outname}"""
