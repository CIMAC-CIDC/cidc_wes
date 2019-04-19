#module: clonality by Pyclone and vcftotsb by vcpileup.seqz.txt.gzftools
#import os
#from string import Template


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


def getNormal_recalsample(wildcards):
    return clonal_trial_runsHelper(wildcards, 0)

def getTumor_recalsample(wildcards):
    return clonal_trial_runsHelper(wildcards, 1)

def clonal_trial_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s.seqz.txt.gz" % (run,run)),
        ls.append("analysis/clonality/%s/%s_sequenza_multibam2seqz.done.txt" % (run,run)),
        ls.append("analysis/clonality/%s/%s.bin50.seqz.txt.gz" % (run,run)),
        ls.append("analysis/clonality/%s/%s.bin50.final.seqz.txt.gz" % (run,run))
        #ls.append("analysis/testclonality/%s/%s.pileup.seqz.txt.gz" % (run,run))
        #ls.append("analysis/clonality/%s/%s.bin50.seqz.txt.gz" % (run,run))
        return ls

rule clonal_trial_all:
    input:
        clonal_trial_targets


rule sequenza_multibam2seqz:
    input:
        tumor_bam=getTumor_recalsample,
        normal_bam=getNormal_recalsample,
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
    #benchmark:
        #"benchmarks/testclonality/{run}/{run}_{run}_preclonality.txt"
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
        headerfile="cidc_wes/header.txt.gz"
    output:
        finalsequenzaoutput="analysis/clonality/{run}/{run}.bin50.final.seqz.txt.gz"
    shell:
        #"zcat {input.headerfile} {input.binned_file_in} |bgzip > {output.finalsequenzaoutput}"
         "cat {input.headerfile} {input.binned_file_in} > {output.finalsequenzaoutput}"
        
##step5:
#to check the order of chromosome
# zcat sample_bin50.out1.20.seqz.txt.gz |cut -f 1 |uniq -c

##step6:
##add header to the file sample_bin50.out1.20.seqz.txt.gz
#chromosome      position        base.ref        depth.normal    depth.tumor     depth.ratio     Af      Bf      zygosity.normal GC.percent      good.reads      AB.normal       AB.tumor        tumor.strand


##step7:
#/home/taing/miniconda3/envs/sequenza/bin/Rscript  /mnt/ssd/wes/cidc_wes/modules/scripts/sequenza.R  /mnt/ssd/wes/analysis/clonality/MDA-Run1-pt1-FF/sample_bin50.out1.20.seqz.txt.gz  /mnt/ssd/wes/analysis/clonality/     sample1.2

##step8:
#run pyclone


