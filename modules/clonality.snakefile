#module: clonality by Pyclone and vcftotsb by vcftools
#import os
#from string import Template

#------------------------------------------------------------------------------
# SET 1 of 2 run helpers: returning sample names
#------------------------------------------------------------------------------
def clonalitypost_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        #tmp.append("analysis/clonality/%s/%s.pyclone.tsv" % (sample_name, sample_name))
        tmp.append(sample_name)
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def clonality_getNormal_samplename(wildcards):
    return clonalitypost_runsHelper(wildcards, 0)

def clonality_getTumor_samplename(wildcards):
    return clonalitypost_runsHelper(wildcards, 1)

#------------------------------------------------------------------------------
# SET 2 of 2 run helpers: returning analysis/align/{sample}/{sample}_recalibrated.bam
#------------------------------------------------------------------------------
def clonal_trial_runsHelper(wildcards, iindex):
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
    #print(tmp)
    return tmp


def clonality_getNormal_recalsample(wildcards):
    return clonal_trial_runsHelper(wildcards, 0)

def clonality_getTumor_recalsample(wildcards):
    return clonal_trial_runsHelper(wildcards, 1)

def clonality_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s.seqz.txt.gz" % (run,run)),
        ls.append("analysis/clonality/%s/%s_sequenza_multibam2seqz.done.txt" % (run,run)),
        ls.append("analysis/clonality/%s/%s.bin50.seqz.txt.gz" % (run,run)),
        ls.append("analysis/clonality/%s/%s.bin50.final.seqz.txt.gz" % (run,run))
        #pyclone output
        ls.append("analysis/clonality/%s/%s_pyclone.tsv" % (run,run)),
        ls.append("analysis/clonality/%s/config.yaml" % run)
    return ls

rule clonality_all:
    input:
        clonality_targets

#------------------------------------------------------------------------------
# START RULES from Clonalitymultithreads.snakefile
#------------------------------------------------------------------------------
rule sequenza_multibam2seqz:
    input:
        tumor_bam=clonality_getTumor_recalsample,
        normal_bam=clonality_getNormal_recalsample,
    output:
        "analysis/clonality/{run}/{run}_sequenza_multibam2seqz.done.txt"
    params:
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
        gc_file=config['gc_file'],
        ref=config['genome_fasta'],
        #chroms= ['chr%' % c for c in list(range(1,22+1))]
        #sequenz_path="/home/aashna/.local/bin", #LEN-add this
        #HARD-CODED- CHANGE
        sequenza_path="/home/taing/miniconda3/envs/sequenza/bin/",
        sequenza_out="analysis/clonality/{run}/{run}.seqz.txt.gz"
    #conda:
    #    "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_sequenza_multibam2seqz.txt"
    threads: 12
    shell:
        "{params.sequenza_path}sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {params.sequenza_out}  --parallel {threads} -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 && touch {output}"

rule mergeChroms:
    input:
        "analysis/clonality/{run}/{run}_sequenza_multibam2seqz.done.txt"
    output:
        "analysis/clonality/{run}/{run}.seqz.txt.gz"
    params:
        chroms="\t".join(["chr%s" % str(i) for i in range(1,23)]),
        name= lambda wildcards: wildcards.run,
	gawk_cmd= "gawk \'{if (NR!=1 && $1 != \"chromosome\") {print $0}}\'" 
    benchmark:
        "benchmarks/clonality/{run}/{run}_mergeChroms.txt"
    run:
        files=" ".join(["analysis/clonality/%s/%s.seqz_%s.txt.gz" % (params.name, params.name, chrom) for chrom in params.chroms.split("\t")])
	#print(params.gawk_cmd)
        #shell("echo zcat {files} |gawk '{if (NR!=1 && $1 != \"chromosome\") {print $0}}' | bgzip > {output}")
	shell("zcat {files} | {params.gawk_cmd} | bgzip > {output}")


rule clonality_sequenza_binning:
    input:
        completeseq="analysis/clonality/{run}/{run}.seqz.txt.gz",
    output:
        sequenzafinal_out="analysis/clonality/{run}/{run}.bin50.seqz.txt.gz"
    params:
        gc_file=config['gc_file'],
        ref=config['genome_fasta'],
        #sequenz_path="/home/aashna/.local/bin", #LEN-
        #HARD-CODED- CHANGE
        sequenza_path="/home/taing/miniconda3/envs/sequenza/bin/"
    #conda:
    #    "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_prestep2clonality.txt"
    shell:
        "{params.sequenza_path}sequenza-utils  seqz_binning --seqz {input.completeseq}  --window 50  -o {output.sequenzafinal_out}"


rule clonality_addheader:
    input:
        binned_file_in="analysis/clonality/{run}/{run}.bin50.seqz.txt.gz",
        headerfile="cidc_wes/static/clonality/header.txt.gz"
    output:
        finalsequenzaoutput="analysis/clonality/{run}/{run}.bin50.final.seqz.txt.gz"
    benchmark:
        "benchmarks/clonality/{run}/{run}_clonality_addheader.txt"
    shell:
        #"zcat {input.headerfile} {input.binned_file_in} |bgzip > {output.finalsequenzaoutput}"
         "cat {input.headerfile} {input.binned_file_in} > {output.finalsequenzaoutput}"
#------------------------------------------------------------------------------
# END RULES from Clonalitymultithreads.snakefile
#------------------------------------------------------------------------------

rule sequenza_fileprep:
    input:
        bin50="analysis/clonality/{run}/{run}.bin50.final.seqz.txt.gz"
    params:
        out_dir="analysis/clonality",
        sample_name=lambda wildcards:wildcards.run,
        #HARD-CODED- CHANGE
        sequenza_env="/home/taing/miniconda3/envs/sequenza/bin/",
    output:
        pyclone_tsv="analysis/clonality/{run}/{run}_pyclone.tsv"
    #conda:
    #    "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}_{run}_postRscript.txt"
    shell:
        #"""Rscript cidc_wes/modules/scripts/sequenza.R  {input}  {params.out_dir}/{wildcards.run}  {params.sample_name}"""
        """{params.sequenza_env}Rscript cidc_wes/modules/scripts/sequenza.R  {input}  {params.out_dir}/{wildcards.run}  {params.sample_name}"""

rule pyclone_finalprocessing:
    input:
        "analysis/clonality/{run}/{run}_pyclone.tsv"
    output:
        "analysis/clonality/{run}/config.yaml"
        #"analysis/clonality/{run}/plots/{run}.pdf",
        #"analysis/clonality/{run}/tables/{run}.tsv",
        #"analysis/clonality/{run}/trace/{run}.txt",
        #"analysis/clonality/{run}/yaml/{run}.yaml",
        #out_dir="analysis/clonality/pyclone"
    params:
        #out_dir=lambda wildcards: "analysis/clonality/%s/" % (wildcards.run, wildcards.run)
        out_dir="analysis/clonality",
        pyclone_env="/home/taing/miniconda3/envs/pyclone/bin/"
    #conda:
    #    "../envs/pyclone.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone.analysis.txt"
    shell:
        """{params.pyclone_env}PyClone run_analysis_pipeline --in_files {input}  --working_dir {params.out_dir}/{wildcards.run} """
