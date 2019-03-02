#module: clonality by Pyclone and vcftools  by vcftools
#import os
#from string import Template


def clonalitypreprocess_runsHelper(wildcards, iindex):
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

def getNormal_sample(wildcards):
    return clonalitypreprocess_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return clonalitypreprocess_runsHelper(wildcards, 1)

def clonalitypreprocess_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s.seqz.txt.gz" % (run,run))
        ls.append("analysis/clonality/%s/%s.bin50.seqz.txt.gz" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone.tsv" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone_analysis/config.yaml" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone_analysis/plots/%s.pdf" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone_analysis/tables/%s.tsv" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone_analysis/trace/%s.txt" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone_analysis/yaml/%s.yaml" % (run,run))
        return ls

rule clonality_all:
    input:
        clonality_targets

rule sequenza_processing:
    input:
        tumor_bam=getTumor_sample,
        normal_bam=getNormal_sample,
    output:
        sequenza_out="analysis/clonality/{run}/{run}.seqz.txt.gz"
    params:
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
        gc_file="/mnt/ssd/wes/cidc_wes/hg38.gc50Base.wig.gz",
        ref=config['genome_fasta'],
        #sequenz_path="/home/aashna/.local/bin", #LEN-add this 
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_preclonality.txt"
    shell:
      "/home/aashna/.local/bin/sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {output.sequenza_out}"
    
rule sequenza_preprocessing:
    input:
        completeseq="analysis/clonality/{run}/{run}.seqz.txt.gz",
    output:
        sequenzafinal_out="analysis/clonality/{run}/{run}.bin50.seqz.txt.gz"
    params:
        #JUST sample names - can also use the helper fns, e.g.                                                                                                                                              
        #normal = lambda wildcards: config['runs'][wildcards.run][0],                                                                                                                                       
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],                                                                                                                                        
        gc_file="/mnt/ssd/wes/cidc_wes/hg38.gc50Base.wig.gz",
        ref=config['genome_fasta'],
        #sequenz_path="/home/aashna/.local/bin", #LEN-add this                                                                                                                                              
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_prestep2clonality.txt"
    shell:
      "/home/aashna/.local/bin/sequenza-utils  seqz_binning --seqz {input.completeseq}  --window 50  -o {output.sequenzafinal_out}"


rule sequenza_fileprep:
  input:
    bin50="analysis/clonality/{run}/{run}.bin50.seq.txt.gz"  #CA209009_run.bin50.seqz.txt.gz
  params:
    out_dir="/analysis/clonality/",
    sample_name= lambda  wildcards: wildcards.run,
output:
      pyclone_tsv="analysis/clonality/{run}/{run}_pyclone.tsv"
  benchmark:
      "benchmarks/clonality/{run}_{run}_postRscript.txt"
  shell:
      """"Rscript  cidc_wes/modules/scripts/sequenza.R  {input}  {params.out_dir}/wildcards.run  {params.samplename}"""



rule pyclone_finalprocessing:
       input:
          intsv="analysis/clonality/{run}/{run}_pyclone.tsv"
       output:
           pyclone_configyaml="analysis/clonality/{run}/{run}_pyclone_analysis/config.yaml",
           pyclone_plots="analysis/clonality/{run}/{run}_pyclone_analysis/plots/{run}.pdf",
           pyclone_tables="analysis/clonality/{run}/{run}_pyclone_analysis/tables/{run}.tsv",
           pyclone_trace="analysis/clonality/{run}/{run}_pyclone_analysis/trace/{run}.txt",
           pyclone_yaml="analysis/clonality/{run}/{run}_pyclone_analysis/yaml/{run}.yaml",
      params:
           out_dir="analysis/clonality/"
     benchmark:
           "benchmarks/clonality/{run}/{run}.pyclone.analysis.txt"
       shell:
           """PyClone run_analysis_pipeline --in_files {input.intsv} --working_dir {params.out_dir}"""


