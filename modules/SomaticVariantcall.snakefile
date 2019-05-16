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
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.vcf.gz" % (run,run))
        ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.output.vcf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.output.vcf.gz" % (run,run))
        #FILTERED VCF
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.filter.vcf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.filter.vcf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.filter.vcf" % (run,run))
        #MAF
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.output.maf" % (run,run))
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
        #Filtered Mutation Signatures
        #ls.append("analysis/somaticVariants/%s/%s_tnsnv.filter.pdf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnhaplotyper2.filter.pdf" % (run,run))
        #ls.append("analysis/somaticVariants/%s/%s_tnscope.filter.pdf" % (run,run))

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
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index}  -i {input.normalbam}  -i  {input.tumorbam}  --algo TNhaplotyper2 --pon {params.tnhaplotyper_pon}  --tumor_sample {params.tumor} --normal_sample {params.normal}   {output.tnhaplotyper2vcf}"""


       

    

    