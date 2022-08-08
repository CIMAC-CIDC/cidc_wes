#module: clonality by Pyclone and vcftotsb by vcftools
#import os
#from string import Template

_clonality_threads=16

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
        #ls.append("analysis/clonality/%s/%s.seqz.txt.gz" % (run,run)),
        #ls.append("analysis/clonality/%s/%s_sequenza_multibam2seqz.done.txt" % (run,run)),
        #ls.append("analysis/clonality/%s/%s.bin50.seqz.txt.gz" % (run,run)),
        #ls.append("analysis/clonality/%s/%s.bin50.final.seqz.txt.gz" % (run,run))
        ls.append("analysis/clonality/%s/%s_genome_view.pdf" % (run,run))
        ls.append("analysis/clonality/%s/%s_segments.txt" % (run,run))
        ls.append("analysis/clonality/%s/%s_sequenza_gainLoss.bed" % (run,run))
	ls.append("analysis/clonality/%s/%s_CP_contours.pdf" % (run,run)) 
	ls.append("analysis/clonality/%s/%s_alternative_solutions.txt" % (run,run))
	ls.append("analysis/clonality/%s/%s_chromosome_view.pdf" % (run,run))
        
        #pyclone output
        #NOTE: _pyclone6.input.tsv should be aggregated across samples for true, multisample clonality analysis
        ls.append("analysis/clonality/%s/%s_pyclone6.input.tsv" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone6.results.tsv" % (run,run))
        ls.append("analysis/clonality/%s/%s_pyclone6.results.summary.tsv" % (run,run))	
        
        #Generate summary json
        ls.append("analysis/report/json/clonality/%s.clonality.json" % run)
    return ls

rule clonality_all:
    input:
        clonality_targets
    benchmark: "benchmarks/clonality/clonality_all.txt"

#------------------------------------------------------------------------------
# START RULES from Clonalitymultithreads.snakefile
#------------------------------------------------------------------------------
rule clonality_sequenza_multibam2seqz:
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
        sequenza_out="%sanalysis/clonality/{run}/{run}.seqz.txt.gz" % config['remote_path'],
        sequenza_bin_path="%s/bin/" % config['sequenza_root'],
    conda:
        "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_sequenza_multibam2seqz.txt"
    threads: 18 #_clonality_threads
    shell:
        "{params.sequenza_bin_path}sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {params.sequenza_out}  --parallel {threads} -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 && touch {output}"

def clonality_mergeChroms_helper(wildcards):
    run = wildcards.run
    chroms=["chr%s" % str(i) for i in range(1,23)]
    files = " ".join(["analysis/clonality/%s/%s.seqz_%s.txt.gz" % (run, run, chrom) for chrom in chroms])
    return files

rule clonality_mergeChroms:
    input:
        "analysis/clonality/{run}/{run}_sequenza_multibam2seqz.done.txt"
    output:
        "analysis/clonality/{run}/{run}.seqz.txt.gz"
    params:
        files = lambda wildcards: clonality_mergeChroms_helper(wildcards),
        name= lambda wildcards: wildcards.run,
	gawk_cmd= "gawk \'{if (NR!=1 && $1 != \"chromosome\") {print $0}}\'" 
    benchmark:
        "benchmarks/clonality/{run}/{run}_mergeChroms.txt"
    shell:
        "zcat {params.files} | {params.gawk_cmd} | bgzip > {output}"


rule clonality_sequenza_binning:
    input:
        completeseq="analysis/clonality/{run}/{run}.seqz.txt.gz",
    output:
        sequenzafinal_out="analysis/clonality/{run}/{run}.bin50.seqz.txt.gz"
    params:
        gc_file=config['gc_file'],
        ref=config['genome_fasta'],
        sequenza_bin_path="%s/bin/" % config['sequenza_root'],
    conda:
        "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_sequenza_binning.txt"
    shell:
        "{params.sequenza_bin_path}sequenza-utils  seqz_binning --seqz {input.completeseq}  --window 50  -o {output.sequenzafinal_out}"


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

# RUN sequenza
rule clonality_sequenza:
    input:
        bin50="analysis/clonality/{run}/{run}.bin50.final.seqz.txt.gz"
    params:
        out_dir="%sanalysis/clonality" % config['remote_path'],
        sample_name=lambda wildcards:wildcards.run,
        sequenza_bin_path="%s/bin/" % config['sequenza_root'],
    output:
        pyclone_tsv="analysis/clonality/{run}/{run}_pyclone.tsv",
        genome_view_plot="analysis/clonality/{run}/{run}_genome_view.pdf",
        segments="analysis/clonality/{run}/{run}_segments.txt",
	CP_contours="analysis/clonality/{run}/{run}_CP_contours.pdf",
	alt_solutions="analysis/clonality/{run}/{run}_alternative_solutions.txt",
	chromosome_view="analysis/clonality/{run}/{run}_chromosome_view.pdf"
	
    conda:
        "../envs/sequenza.yml"
    benchmark:
        "benchmarks/clonality/{run}/{run}_sequenza.txt"
    shell:
        """{params.sequenza_bin_path}Rscript cidc_wes/modules/scripts/sequenza.R  {input}  {params.out_dir}/{wildcards.run}  {params.sample_name}"""

rule clonality_sequenza2pyclone6:
    input:
        "analysis/clonality/{run}/{run}_pyclone.tsv",
    params:
        run = lambda wildcards: wildcards.run,
    benchmark:
        "benchmarks/clonality/{run}/{run}_sequenza2pyclone6.txt"
    output:
        "analysis/clonality/{run}/{run}_pyclone6.input.tsv",
    shell:
        "./cidc_wes/modules/scripts/sequenza2pyclone6.py -f {input} -n {params.run} -o {output}"

rule clonality_pyclone6_fit:
    input:
        "analysis/clonality/{run}/{run}_pyclone6.input.tsv"
    params:
        num_clusters = 40, #between 10-40;
        density = "beta-binomial", #alt: binomial
        #num-grid-points = 100, #default: 100, use higher for deeply sequenced data
        num_restarts = 10, #between 10-100, higher val means longer run-time
    benchmark:
        "benchmarks/clonality/{run}/{run}_pyclone6_fit.txt"
    output:
        "analysis/clonality/{run}/{run}_pyclone6.h5",
    shell:
        "pyclone-vi fit -i {input} -o {output} -c {params.num_clusters} -d {params.density} -r {params.num_restarts}"

rule clonality_pyclone6_writeOut:
    input:
        "analysis/clonality/{run}/{run}_pyclone6.h5"
    benchmark:
        "benchmarks/clonality/{run}/{run}_pyclone6_writeOut.txt"
    output:
        "analysis/clonality/{run}/{run}_pyclone6.results.tsv"
    shell:
        "pyclone-vi write-results-file -i {input} -o {output}"

rule clonality_pyclone6_summarizeResults:
    input:
        "analysis/clonality/{run}/{run}_pyclone6.results.tsv"
    benchmark:
        "benchmarks/clonality/{run}/{run}_pyclone6_summarizeResults.txt"
    output:
        "analysis/clonality/{run}/{run}_pyclone6.results.summary.tsv"
    shell:
        #remove hdr
        #"cut -f 3,4 {input} | tail -n +2 | uniq > {output}"
        #keep hdr
        "cut -f 3-5 {input} | uniq > {output}"

rule clonality_json:
    """jsonify the tumor cloanlity
    """
    input:
        #NOTE: 2 output files _pyclone6.results.tsv and .summary.tsv
        #using .summary.tsv for json summary, but should ingest .results.tsv
        results="analysis/clonality/{run}/{run}_pyclone6.results.tsv",
        summary="analysis/clonality/{run}/{run}_pyclone6.results.summary.tsv",
        pyclone6_input="analysis/clonality/{run}/{run}_pyclone6.input.tsv",
    output:
        "analysis/report/json/clonality/{run}.clonality.json"
    params:
        run = lambda wildcards: wildcards.run
    group: "clonality"
    benchmark:
        "benchmarks/clonality/{run}.clonality_json.txt"
    shell:
        "cidc_wes/modules/scripts/json_clonality.py -r {params.run} -i {input.pyclone6_input} -j {input.results} -k {input.summary} -o {output}"

rule clonality_callGainLoss:
    """use hard-cutoffs to call regions of GAIN/LOSS"""
    input:
        "analysis/clonality/{run}/{run}_segments.txt"
    output:
        #NOTE: changing from {run}-{tmr} to {run}-{run} to be more consistent
        "analysis/clonality/{run}/{run}_sequenza_gainLoss.bed"
    group: "clonality"
    log: "analysis/logs/clonality/{run}/{run}.callGainLoss.log"
    benchmark:
        "benchmarks/clonality/{run}/{run}.callGainLoss.txt"
    shell:
        "./cidc_wes/modules/scripts/copynumber_callGainLoss.py -f {input} -o {output}"

# OBSOLETE
# rule pyclone_build_mutation_file:
#     input:
#         "analysis/clonality/{run}/{run}_pyclone.tsv"
#     params:
#         pyclone_bin_path="%s/bin/" % config['pyclone_root'],
#     output:
#         "analysis/clonality/{run}/{run}_pyclone.yaml"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone_build_mutation_file.txt"
#     shell:
#         "{params.pyclone_bin_path}PyClone build_mutations_file --in_file {input} --out_file {output}"

# rule pyclone_generate_config:
#     input:
#         "analysis/clonality/{run}/{run}_pyclone.yaml"
#     output:
#         "analysis/clonality/{run}/pyclone.config.yaml"
#     params:
#         outdir = lambda wildcards: "%sanalysis/clonality/%s" % (config['remote_path'],wildcards.run),
#         template = "cidc_wes/static/clonality/pyclone.config.yaml"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone.generate.config.txt"
#     shell:
#         "cidc_wes/modules/scripts/clonality_pycloneConfig.py -t {params.template} -m {input} -o {params.outdir} > {output}"

# rule pyclone_run_analysis:
#     input:
#         "analysis/clonality/{run}/pyclone.config.yaml"
#     params:
#         pyclone_bin_path="%s/bin/" % config['pyclone_root'],
#     output:
#         #LEN: I don't know what to key in on as output
#         "analysis/clonality/{run}/trace/alpha.tsv.bz2"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone.run_analysis.txt"
#     shell:
#         "{params.pyclone_bin_path}PyClone run_analysis --config_file {input}"

# rule pyclone_build_table:
#     input:
#         conf="analysis/clonality/{run}/pyclone.config.yaml",
#         other="analysis/clonality/{run}/trace/alpha.tsv.bz2",
#     params:
#         pyclone_bin_path="%s/bin/" % config['pyclone_root'],
#     output:
#         "analysis/clonality/{run}/{run}_table.tsv"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone.build_table.txt"
#     shell:
#         "{params.pyclone_bin_path}PyClone build_table --config_file {input.conf} --out_file {output} --table_type cluster --max_clusters 100"

# rule pyclone_density_plot:
#     input:
#         conf="analysis/clonality/{run}/pyclone.config.yaml",
#         other="analysis/clonality/{run}/trace/alpha.tsv.bz2",
#     params:
#         pyclone_bin_path="%s/bin/" % config['pyclone_root'],
#     output:
#         "analysis/clonality/{run}/{run}_plot.density.pdf"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone_density_plot.txt"
#     shell:
#         "{params.pyclone_bin_path}PyClone plot_clusters --config_file {input.conf} --plot_file {output} --plot_type density --max_clusters 100"

# rule pyclone_scatter_plot:
#     input:
#         conf="analysis/clonality/{run}/pyclone.config.yaml",
#         other="analysis/clonality/{run}/trace/alpha.tsv.bz2",
#     params:
#         pyclone_bin_path="%s/bin/" % config['pyclone_root'],
#     output:
#         "analysis/clonality/{run}/{run}_plot.scatter.pdf"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone_scatter_plot.txt"
#     shell:
#         "{params.pyclone_bin_path}PyClone plot_clusters --config_file {input.conf} --plot_file {output} --plot_type scatter --max_clusters 100"

    
# rule pyclone_parallelCoordinates_plot:
#     input:
#         conf="analysis/clonality/{run}/pyclone.config.yaml",
#         other="analysis/clonality/{run}/trace/alpha.tsv.bz2",
#     params:
#         pyclone_bin_path="%s/bin/" % config['pyclone_root'],
#     output:
#         "analysis/clonality/{run}/{run}_plot.coordinates.pdf"
#     conda:
#         "../envs/pyclone.yml"
#     benchmark:
#         "benchmarks/clonality/{run}/{run}.pyclone_coordinates_plot.txt"
#     shell:
#         "{params.pyclone_bin_path}PyClone plot_clusters --config_file {input.conf} --plot_file {output} --plot_type parallel_coordinates --max_clusters 100"

