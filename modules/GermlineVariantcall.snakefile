# #module: germline calls by Sentieon,added all samples
# #import os
# #from string import Template
# _germlinecalls_threads=32
# _vcf2maf_threads=8


# #NOT used 
# # #NOTE: germline_runsHelper, getNormal_sample, and germ_getTumor are NOT
# # #called by any one!
# # def germline_runsHelper(wildcards, iindex):
# #     """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
# #     returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
# #     tmp = []
# #     r = config['runs'][wildcards.run]
# #     #print(r)

# #     #check that we have a valid pair
# #     if len(r) >=2:
# #         sample_name = r[iindex]
# #         tmp.append("analysis/align/%s/%s_recalibrated.bam" % (sample_name, sample_name))
# #     else:
# #         #NOTE: I can't figure out a proper kill command so I'll do this
# #         tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
# #     #print(tmp)
# #     return tmp


# # def germ_getNormal(wildcards):
# #     return germline_runsHelper(wildcards, 0)

# # def germ_getTumor(wildcards):
# #     return germline_runsHelper(wildcards, 1)

# def germlinecalls_targets(wildcards):
#     """Generates the targets for this module"""
#     ls = []
#     for sample in config['samples']:
#         #Consolidate these with an inner-for-loop?
#         ls.append("analysis/germlineVariants/%s/%s_dnascope.output.vcf.gz" % (sample,sample))
#         ls.append("analysis/germlineVariants/%s/%s_haplotyper.output.vcf.gz" % (sample,sample))
# 	#FILTERED VCF
#         ls.append("analysis/germlineVariants/%s/%s_dnascope.output.filter.vcf" % (sample,sample))
#         ls.append("analysis/germlineVariants/%s/%s_haplotyper.output.filter.vcf" % (sample,sample))
#         #MAF
#         ls.append("analysis/germlineVariants/%s/%s_dnascope.output.maf" % (sample,sample))
#         ls.append("analysis/germlineVariants/%s/%s_haplotyper.output.maf" % (sample,sample))
#         #READ DEPTH/COVERAGE FILTER: 10x,30x
#         for frac in [10, 30]:
#             ls.append("analysis/germlineVariants/%s/%s_dnascope.coverage.%s.vcf" % (sample,sample, str(frac)))
#             ls.append("analysis/germlineVariants/%s/%s_haplotyper.coverage.%s.vcf" % (sample,sample, str(frac)))
#             #VCF-COMPARISON
#             ls.append("analysis/germlineVariants/%s/%s_comparedsamples.diff.discordance_matrix" % (sample,sample))
#             return ls

# rule germlinecalls_all:
#     input:
#         germlinecalls_targets
	
# rule germline_calling_DNAscope:
#     input:
#         in_recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam"
#     output:
#         dnascopevcf="analysis/germlineVariants/{sample}/{sample}_dnascope.output.vcf.gz"
#     params:
#         index=config['genome_fasta'],
#         sentieon_path=config['sentieon_path'],
#         dbsnp= config['dbsnp'],
#         #JUST sample names - can also use the helper fns, e.g.
#         #normal = lambda wildcards: config['runs'][wildcards.run][0],
#         #tumor = lambda wildcards: config['runs'][wildcards.run][1],
#     threads:_germlinecalls_threads
#     benchmark:
#         "benchmarks/germlineVariantscall/{sample}/{sample}.germline_calling_DNAscope.txt"
#     shell:
#         """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.in_recalibratedbam} --algo DNAscope --dbsnp {params.dbsnp}  --emit_conf=30 --call_conf=30 {output.dnascopevcf}"""


# rule germline_calling_haplotyper:
#     input:
#         in_recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam"
#     output:
#         haplotypervcf="analysis/germlineVariants/{sample}/{sample}_haplotyper.output.vcf.gz"
#     params:
#         index=config['genome_fasta'],
#         sentieon_path=config['sentieon_path'],
#         dbsnp= config['dbsnp'],
#         #JUST sample names - can also use the helper fns, e.g.
#         #normal = lambda wildcards: germ_getNormal(wildcards)
#         #normal = lambda wildcards: config['runs'][wildcards.run][0],
#         #tumor = lambda wildcards: config['runs'][wildcards.run][1],
#     threads:_germlinecalls_threads
#     benchmark:
#         "benchmarks/germlineVariantscall/{sample}/{sample}.germline_calling_haplotyper.txt"
#     shell:
#         """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.in_recalibratedbam} --algo Haplotyper  --dbsnp {params.dbsnp}  --emit_conf=30 --call_conf=30 {output.haplotypervcf}"""


# rule germline_vcftoolsfilter:
#     """General rule to filter the two different types of vcf.gz files"""
#     input:
#         vcffiles="analysis/germlineVariants/{run}/{run}_{caller}.output.vcf.gz"
#     output:
#         filteredvcf="analysis/germlineVariants/{run}/{run}_{caller}.output.filter.vcf"
#     params:
#         index=config['genome_fasta'],
#     benchmark:
#         "benchmarks/germlineVariantscall/{run}/{run}.{caller}_vcftoolsfilter.txt"
#     shell:
#         """vcftools --gzvcf {input.vcffiles} --remove-filtered-all --recode --stdout > {output.filteredvcf}"""


# rule germline_gunzip_vcf:
#     """General rule to gunzip the two  types of vcf.gz files-DNAscope and haplotyper"""
#     input:
#         "analysis/germlineVariants/{run}/{run}_{caller}.output.vcf.gz"
#     output:
#         #Should we make this a temp?
#         "analysis/germlineVariants/{run}/{run}_{caller}.output.vcf"
#     benchmark:
#         "benchmarks/germlineVariantscall/{run}/{run}.{caller}_gunzip_vcf.txt"
#     shell:
#         #NOTE: we want to keep the original .gz vcf file
#         "gunzip -k {input}"

# rule germline_vcf2maf:
#     """General rule to convert the different vcf files into maf"""
#     input:
#         "analysis/germlineVariants/{run}/{run}_{caller}.output.vcf"
#     output:
#         "analysis/germlineVariants/{run}/{run}_{caller}.output.maf"
#     threads: _vcf2maf_threads,
#     params:
#         index=config['genome_fasta'],
#         vep_path="%s/bin" % config['wes_root'],
#         vep_data=config['vep_data'],
#         vep_assembly=config['vep_assembly'],
#     benchmark:
#         "benchmarks/germlineVariantscall/{run}/{run}.{caller}_vcf2maf.txt"
#     shell:
#         """vcf2maf.pl --input-vcf {input} --output-maf {output} --ref-fasta {params.index} --vep-path {params.vep_path} --vep-data {params.vep_data} --ncbi-build {params.vep_assembly}"""

# rule coverage_filter_DNAscope:
#     input:
#         "analysis/germlineVariants/{sample}/{sample}_dnascope.output.vcf.gz"
#     params:
#         threshold=lambda wildcards: wildcards.frac,
#         field="DP" #NOTE this is the particular field germline  vcf files
#     output:
#         #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
#         #with vcftoolsfilter; {frac} is int
#         "analysis/germlineVariants/{sample}/{sample}_dnascope.coverage.{frac,\d+}.vcf"
#     benchmark:
#         "benchmarks/germlineVariantscall/{sample}/{sample}.coverage_filter_dnascope.txt"
#     shell:
#         "cidc_wes/modules/scripts/vcf_filterByReadDepth.py -v {input} -t {params.threshold} -f {params.field} -o {output}"

# rule coverage_filter_germlinehaplotyper:
#     input:
#         "analysis/germlineVariants/{sample}/{sample}_haplotyper.output.vcf.gz"
#     params:
#         threshold=lambda wildcards: wildcards.frac,
#         field="DP" #NOTE this is the particular field germline haplotyper vcf files
#     output:
#         #NOTE: need to add regular-expression for {frac} b/c it's ambiguous
#         #with vcftoolsfilter; {frac} is int
#         "analysis/germlineVariants/{sample}/{sample}_haplotyper.coverage.{frac,\d+}.vcf"
#     benchmark:
#         "benchmarks/germlineVariantscall/{sample}/{sample}.coverage_filter_haplotyper.txt"
#     shell:
#         "cidc_wes/modules/scripts/vcf_filterByReadDepth.py -v {input} -t {params.threshold} -f {params.field} -o {output}"

# rule filterOutRandomContigs:
#     """NOTE: for vcfintersect_bedtools to work properly, the files must 
#     agree on the same chromosomes.  
#     For example, if chr1_KI270709v1_random is in one vcf file, it must 
#     also be in the other.  
#     OTHERWISE --diff-discordance-matrix will give the following error-
#     Found chr1_KI270709v1_random in file 1 and chr1_KI270711v1_random in file 2

#     NOTE: it recommends using --not-chr, but that requires listing out
#     the chromosomes to EXCLUDE and in hg38 there are many randome ones!

#     It's easier to just filter out everything except chr1-23,X,Y,M"""
#     input:
#       "analysis/germlineVariants/{run}/{run}_{caller}.output.vcf"
#     output:
#       "analysis/germlineVariants/{run}/{run}_{caller}.canonical.vcf"
#     params:
#     #got this grep cmd from here-
#     #ref: https://www.biostars.org/p/201603/
#       grep_cmd = "\'^#\|^#CHROM\|^chr[1-23,X,Y,M]\'" #HARD-code chr1-23,X,T,M
#     shell:
#       "grep -w {params.grep_cmd} {input} > {output}"


# rule vcfintersect_bedtools:
#     input:
#         #dnascopevcf="analysis/germlineVariants/{sample}/{sample}_dnascope.output.vcf",
#         dnascopevcf="analysis/germlineVariants/{sample}/{sample}_dnascope.canonical.vcf",
#         #haplotypervcf="analysis/germlineVariants/{sample}/{sample}_haplotyper.output.vcf"
#         haplotypervcf="analysis/germlineVariants/{sample}/{sample}_haplotyper.canonical.vcf"
#     output:
#         comparedfiles="analysis/germlineVariants/{sample}/{sample}_comparedsamples.diff.discordance_matrix"
#     params:
#         outfile=lambda wildcards: "analysis/germlineVariants/%s/%s_comparedsamples" % (wildcards.sample, wildcards.sample)
#     benchmark:
#         "benchmarks/germlineVariantscall/{sample}/{sample}.vcfintersect_comparedsamples.txt"
#     shell:
#         #NOTE to aashna: the error with using --vcf instead of --diff for the
#         #second input file
#         #"""vcftools --diff-discordance-matrix  --vcf {input.dnascopevcf}  --vcf {input.haplotypervcf}  --out {output.comparedfiles}"""
#         """vcftools --diff-discordance-matrix  --vcf {input.dnascopevcf}  --diff {input.haplotypervcf}  --out {params.outfile}"""


# module: Generate SNP92 to identification of germline variant.
# paper_name:SNPs for a universal individual identification panel

def germlinecalls_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        #Consolidate these with an inner-for-loop?
        ls.append("analysis/germlineVariants/%s/%s_variant.vcf" % (sample,sample))
        ls.append("analysis/germlineVariants/%s/%s_SNP92.recode.vcf" % (sample,sample))

rule germlinecalls_all:
    input:
        germlinecalls_targets

rule germlinecalls_bamtovcf:
    input:
        input_sortbamfile = "analysis/align/{sample}/{sample}_sorted.bam" 
    output:
        output_vcf="analysis/germlineVariants/{sample}/{sample}_variant.vcf"
    params:
        index=config['genome_fasta']
        positions_bamtovcf=['positions_bamtovcf']
    shell:
        """samtools mpileup -g -Q 0 -f index {input.input_sortbamfile} --positions {params.positions_bamtovcf} | bcftools view > {output.output_vcf}"""

rule germlinecalls_snp92:
    input:
        input_vcf="analysis/germlineVariants/{sample}/{sample}_variant.vcf"
    output:
        output_SNP92="analysis/germlineVariants/{sample}/{sample}_SNP92"
    params:
        positons_SNP92=config['positions_SNP92']
    shell:
        """vcftools --vcf {input.input_vcf} --positions {params.positons_SNP92} --recode --out {output.output_SNP92}"""
