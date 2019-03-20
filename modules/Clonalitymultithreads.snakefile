#module: clonality by Pyclone and vcftotsb by vcpileup.seqz.txt.gzftools
#import os
#from string import Template





###step1:
#add this to align module
#samtools sort *_recalibrated.bam > *sort_recalibrated.bam
#samtools index *sort_recalibrated.bam



###step2:

def clonal_trial_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/align/%s/%s.sort_recalibrated.bam" % (sample_name, sample_name))
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
        ls.append("analysis/testclonality1/%s/%s.seqz.txt.gz" % (run,run))
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
        sequenza_out="analysis/testclonality1/{run}/{run}.seqz.txt.gz"
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
        #sequenza_out="analysis/testclonality1"
    #conda:
    #    "../envs/sequenza.yml"
    #benchmark:
        #"benchmarks/testclonality/{run}/{run}_{run}_preclonality.txt"
    shell:
        "{params.sequenza_path}sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {output.sequenza_out}  --parallel 12 -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 ch16 chr17 chr18 chr19 chr20 chr21 chr22 "



 ##step3:
#" zcat MDA-Run1-pt1-FF.seqz_chr1.txt.gz    MDA-Run1-pt1-FF.seqz_chr2.txt.gz MDA-Run1-pt1-FF.seqz_chr3.txt.gz  MDA-Run1-pt1-FF.seqz_chr4.txt.gz  MDA-Run1-pt1-FF.seqz_chr5.txt.gz MDA-Run1-pt1-FF.seqz_chr6.txt.gz MDA-Run1-pt1-FF.seqz_chr7.txt.gz MDA-Run1-pt1-FF.seqz_chr8.txt.gz  MDA-Run1-pt1-FF.seqz_chr9.txt.gz  MDA-Run1-pt1-FF.seqz_chr10.txt.gz MDA-Run1-pt1-FF.seqz_chr11.txt.gz MDA-Run1-pt1-FF.seqz_chr12.txt.gz MDA-Run1-pt1-FF.seqz_chr13.txt.gz MDA-Run1-pt1-FF.seqz_chr14.txt.gz MDA-Run1-pt1-FF.seqz_ch16.txt.gz MDA-Run1-pt1-FF.seqz_chr17.txt.gz MDA-Run1-pt1-FF.seqz_chr18.txt.gz MDA-Run1-pt1-FF.seqz_chr19.txt.gz MDA-Run1-pt1-FF.seqz_chr20.txt.gz MDA-Run1-pt1-FF.seqz_chr21.txt.gz MDA-Run1-pt1-FF.seqz_chr22.txt.gz  |gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip > sample.out.seqz.txt.gz"

#To be noted - chromosome 16 file name has a bug in the software . File name comes as "ch16"


###step4:
#"/home/taing/miniconda3/envs/sequenza/bin/sequenza-utils seqz_binning --window 50  -s  sample.seqz.out1.gz | gzip > sample_bin50.out1.20.seqz.txt.gz"


##step5:
#to check the order of chromosome
# zcat sample_bin50.out1.20.seqz.txt.gz |cut -f 1 |uniq -c

##step6:
##add header to the file sample_bin50.out1.20.seqz.txt.gz
#chromosome      position        base.ref        depth.normal    depth.tumor     depth.ratio     Af      Bf      zygosity.normal GC.percent      good.reads      AB.normal       AB.tumor        tumor.strand

##step7:
#/home/taing/miniconda3/envs/sequenza/bin/Rscript  /mnt/ssd/wes/cidc_wes/modules/scripts/sequenza.R  /mnt/ssd/wes/analysis/testclonality1/MDA-Run1-pt1-FF/sample_bin50.out1.20.seqz.txt.gz  /mnt/ssd/wes/analysis/testclonality1/     sample1.2

##step8:
#run pyclone

