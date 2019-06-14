#MODULE: Purity by Facets
#import os
_puritycalls_threads=16

def puritybam_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/align/%s/%s_recalibrated.bam" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    return tmp


def puritybai_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/align/%s/%s_recalibrated.bam.bai" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    return tmp


def purity_getNormal_bam_sample(wildcards):
    return puritybam_runsHelper(wildcards, 0)

def purity_getTumor_bam_sample(wildcards):
    return puritybam_runsHelper(wildcards, 1)
    
def purity_getNormal_bai_sample(wildcards):
    return puritybai_runsHelper(wildcards, 0)

def purity_getTumor_bai_sample(wildcards):
    return puritybai_runsHelper(wildcards, 1)

    
def purity_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"]:
    	ls.append("analysis/purity/%s/%s_purity_results.txt" % (run,run))

        #purity post process
        ls.append("analysis/purity/%s/%s_purity_postprocessed_results.txt" % (run,run))
        ls.append("analysis/purity/%s/%s.cncf" % (run,run))
        ls.append("analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run))
        ls.append("analysis/purity/%s/%s.iterpurityvalues.txt" % (run,run))

    return ls

rule purity_all:
    input:
        purity_targets

rule Puritycalls_Facets:
    """Get the  recalibrated bam files from  mapped reads"""
    input:
         NormalBam=purity_getNormal_bam_sample,
         TumorBam=purity_getTumor_bam_sample,
         Normalbai=purity_getNormal_bai_sample,
         Tumorbai=purity_getTumor_bai_sample
    output:
         puritycalls="analysis/purity/{run}/{run}_purity_results.txt",
    message:
        "FACETS: purity calculations for bam file"
    group: "purity"
    params:
        index=config['genome_fasta'],
        index1=config['facets_vcftar']
    threads: _puritycalls_threads
    conda:
        "../envs/purity.yml"
    benchmark:
        "benchmarks/puritycalls/{run}/{run}.purityresults.txt"
    shell:
        "snp-pileup -q15 -Q20  {params.index1} {output} {input.NormalBam} {input.TumorBam}"

rule purityprocessing_filter:
    input:
        "analysis/purity/{run}/{run}_purity_results.txt"
    params:
        #no params
    output:
        "analysis/purity/{run}/{run}_purity_postprocessed_results.txt"
    group: "purity"
    benchmark:
        "benchmarks/puritycalls/{run}/{run}.postprocessed.purity.results.txt"
    shell:
        #"cidc_wes/modules/scripts/vcf_filterByReadDepth.py -v {input} -t {params.threshold} -f {params.field} -o {output}"
        """cat {input} | sed 's/chr//g'  > {output} """

rule purityplots_postprocessing:
    input:
        "analysis/purity/{run}/{run}_purity_postprocessed_results.txt"
    params:
        name=lambda wildcards: wildcards.run,
        output_dir=lambda wildcards: "%sanalysis/purity/%s/" % (config['remote_path'], wildcards.run)
    output:
        cncf="analysis/purity/{run}/{run}.cncf",
        opt="analysis/purity/{run}/{run}.optimalpurityvalue.txt",
        iter="analysis/purity/{run}/{run}.iterpurityvalues.txt",
    group: "purity"
    conda:
        "../envs/purity.yml"
    benchmark:
        "benchmarks/puritycalls/{run}/{run}.purity.postprocessingplots.txt"
    shell:
        "Rscript --vanilla cidc_wes/modules/scripts/facets_plots.R {input} {params.output_dir} {params.name}"
