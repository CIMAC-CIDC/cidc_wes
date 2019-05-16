#module: Somatic Variant calls by Sentieon
#import os
#from string import Template

_somaticcall_threads=8
#_vcf2maf_threads=4

#NOTE: somatic_runsHelper, getNormal_sample, and getTumor_sample are NOT
#called by any one!
def somatic_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append(sample_name)
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def somatic_getNormal(wildcards):
    return somatic_runsHelper(wildcards, 0)

def somatic_getTumor(wildcards):
    return somatic_runsHelper(wildcards, 1)

def somaticall_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #Consolidate these with an inner-for-loop?
        #ls.append("analysis/somaticVariants/%s/%s_call.output.stats" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.vcf.gz" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.output.vcf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.output.vcf.gz" % (run,run))
        #FILTERED VCF
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.filter.vcf" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.filter.vcf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.filter.vcf" % (run,run))
        #MAF
        ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.maf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.output.maf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.output.maf" % (run,run))
        #Filtered MAF
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.filter.maf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.filter.maf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.filter.maf" % (run,run))
        #Mutation Signatures
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.pdf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.output.pdf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.output.pdf" % (run,run))
        #EXON mutations- should this be on full or filtered?
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.exon.maf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.output.exon.maf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.output.exon.maf" % (run,run))
        #alleleFrac cutoffs - should this be on full or filtered?
        #for frac in [0.05,0.1,0.2,0.3,0.4,0.5]:
          #ls.append("analysis/somaticVariants/%s/%s_tnscope.output.%s.vcf" % (run,run, str(frac)))
          #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper.output.%s.vcf" % (run,run, str(frac)))

        #read depth/coverage filter: 10x, 20x, 50x - should this be on full or filtered?
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.coverage.%s.vcf" % (run,run, str(frac)))

        #Mutation load
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.mutationload.txt" % (run,run))
        return ls

rule somaticcalls_all:
     input:	
	somaticall_targets

rule somatic_calling_TNhaplotyper2:
     input:
	tumorbam="analysis/align/R1-3-F/R1-3-F_recalibrated.bam",
	tumorbai="analysis/align/R1-3-F/R1-3-F_recalibrated.bam.bai",
	normalbam="analysis/align/R1-3-N/R1-3-N_recalibrated.bam",
	normalbai="analysis/align/R1-3-N/R1-3-N_recalibrated.bam.bai"
     output:
	tnhaplotyper2vcf="analysis/somaticVariants/{run}/{run}_tnhaplotyper2.output.vcf"
     params:
	index=config['genome_fasta'],
	sentieon_path=config['sentieon_path'],
	tnhaplotyper_pon= config['pons_haplotyper'],
	#JUST sample names - can also use the helper fns, e.g.
	#normal = lambda wildcards: getNormal_sample(wildcards)
	normal = lambda wildcards: config['runs'][wildcards.run][0],
	tumor = lambda wildcards: config['runs'][wildcards.run][1],
     threads:_somaticcall_threads
     benchmark:
	"benchmarks/somaticvariantcall/{run}/{run}.somatic_calling_TNhaplotyper2.txt"
     shell:
	"""{params.sentieon_path}/sentieon driver -t {threads} -r {params.index}  -i {input.normalbam}  -i  {input.tumorbam}  --algo TNhaplotyper2 --pon {params.tnhaplotyper_pon}  --tumor_    sample {params.tumor} --normal_sample {params.normal}   {output.tnhaplotyper2vcf}"""


rule somatic_calling_TNsnv:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        statscall="analysis/somaticVariants/{run}/{run}_call.output.stats",
        tnsnvvcf="analysis/somaticVariants/{run}/{run}_tnsnv.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        #JUST sample names - can also use the helper fns, e.g.
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:96
    benchmark:
        "benchmarks/somaticvariantcall/{run}/{run}.somatic_calling_TNsnv.txt"
    shell:
        #"""{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.corealig#nedbam} --algo TNsnv --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.#dbsnp} --call_stats_out {output.statscall} --min_tumor_allele_frac 0.05 {output.tnsnvvcf}"""
        #REMOVING min_tumor_allele_frac param
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo TNsnv --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} --call_stats_out {output.statscall} {output.tnsnvvcf}"""


rule vcftoolsfilter:
     """General rule to filter the three different types of vcf.gz files"""
     input:
	"analysis/somaticVariants/{run}/{run}_{caller}.output.vcf.gz"
     output:
        "analysis/somaticVariants/{run}/{run}_{caller}.filter.vcf"
     params:
        index=config['genome_fasta'],
     benchmark:
        "benchmarks/somaticvariantcall/{run}/{run}.{caller}_vcftoolsfilter.txt"
     shell:
       """vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""


rule gunzip_vcf:
     """General rule to gunzip the three types of vcf.gz files-
        tnscope_, tnsnv, and tnhaplotyper"""
      input:
	"analysis/somaticVariants/{run}/{run}_{caller}.output.vcf.gz"
      output:
	#Should we make this a temp?
        "analysis/somaticVariants/{run}/{run}_{caller}.output.vcf"
      benchmark:
         "benchmarks/somaticvariantcall/{run}/{run}.{caller}_gunzip_vcf.txt"
      shell:
         #NOTE: we want to keep the original .gz vcf file
         "gunzip -k {input}"

rule vcfVEP:
      """Rule to annotate vcf files with vep"""
     input:
         "analysis/somaticVariants/{run}/{run}_{caller}.{type}.vcf"
       output:
         "analysis/somaticVariants/{run}/{run}_{caller}.{type}.vep.vcf"
      params:
          vep_data=config['vep_data'],
          vep_synonyms=config['vep_synonyms'],
      benchmark:
           "benchmarks/somaticvariantcall/{run}/{run}.{caller}.{type}_vcfVEP.txt"
      shell:
           "vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs"


rule vcf2maf:
       """General rule to convert the different vcf files into maf"""
       input:
	 vcf="analysis/somaticVariants/{run}/{run}_{caller}.{type}.vcf",
	 vep="analysis/somaticVariants/{run}/{run}_{caller}.{type}.vep.vcf",
       output:
         "analysis/somaticVariants/{run}/{run}_{caller}.{type}.maf"
       params:
          vep_index=config['vep_fasta'],
          vep_custom_enst= config['vep_custom_enst'],
          vep_assembly=config['vep_assembly'],
          vep_filter= config['vep_filter'],
          buffer_size=config['vcf2maf_bufferSize'],
	  tumor= lambda wildcards: somatic_getTumor(wildcards),
          normal= lambda wildcards: somatic_getNormal(wildcards),
       benchmark:
          "benchmarks/somaticvariantcall/{run}/{run}.{caller}.{type}_vcf2maf.txt"
       log:
          "analysis/log/somaticvariantcall/{run}/{run}.{caller}.{type}_vcf2maf.log.txt"
       shell:
           """vcf2maf.pl --input-vcf {input} --output-maf {output} --custom-enst {params.vep_custom_enst} --ref-fasta {params.vep_index} --tumor-id {params.tumor} --normal-id {params.normal} --ncbi-build {params.vep_assembly} --filter-vcf {params.vep_filter} --buffer-size {params.buffer_size} 2> {log}"""

       
