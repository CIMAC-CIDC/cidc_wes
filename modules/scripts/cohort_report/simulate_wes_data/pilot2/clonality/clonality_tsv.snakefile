def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s_pyclone.tsv" % (run,run))
        ls.append("analysis/clonality/%s/trace/alpha.tsv.bz2" % run)
        #tmr = config[run]['tumor']
        #nrm = config[run]['normal']
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (tmr,tmr))
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (nrm,nrm))


    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

configfile: "config.cohort_json.small.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

def sequenza_multibam2seqzInputfn(wildcards):
    run = wildcards.run
    tmr = config[run]['tumor']
    nrm = config[run]['normal']
    tmp = {}
    tmp['tumor_bam'] = "recal/%s_recalibrated.bam" % tmr
    tmp['normal_bam'] = "recal/%s_recalibrated.bam" % nrm
    #print(tmp)
    return tmp
    
rule sequenza_multibam2seqz:
    input:
        unpack(sequenza_multibam2seqzInputfn)
    output:
        "analysis/clonality/{run}/{run}_sequenza_multibam2seqz.done.txt"
    params:
        #JUST sample names - can also use the helper fns, e.g.
        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
        #gc_file=config['gc_file'],
        #ref=config['genome_fasta'],
        gc_file="./ref_files/hg38/clonality/hg38.gc50Base.txt.gz",
        ref="./ref_files/hg38/bwa_gdc/GRCh38.d1.vd1.CIDC.fa",
        #chroms= ['chr%' % c for c in list(range(1,22+1))]
        sequenza_out="analysis/clonality/{run}/{run}.seqz.txt.gz",
        sequenza_bin_path="~/miniconda3/envs/sequenza/bin/",
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_sequenza_multibam2seqz.txt"
    threads: 16
    shell:
        "{params.sequenza_bin_path}sequenza-utils  bam2seqz -n {input.normal_bam}  -t {input.tumor_bam}  --fasta {params.ref} -gc {params.gc_file} -o {params.sequenza_out}  --parallel {threads} -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 && touch {output}"

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
        gc_file="./ref_files/hg38/clonality/hg38.gc50Base.txt.gz",
        ref="./ref_files/hg38/bwa_gdc/GRCh38.d1.vd1.CIDC.fa",
        sequenza_bin_path="~/miniconda3/envs/sequenza/bin/",
    benchmark:
        "benchmarks/clonality/{run}/{run}_{run}_prestep2clonality.txt"
    shell:
        "{params.sequenza_bin_path}sequenza-utils  seqz_binning --seqz {input.completeseq}  --window 50  -o {output.sequenzafinal_out}"

rule clonality_addheader:
    input:
        binned_file_in="analysis/clonality/{run}/{run}.bin50.seqz.txt.gz",
        #hard-coding
        headerfile="/mnt/ssd/mda-r1-pt1_report/cidc_wes/static/clonality/header.txt.gz"
    output:
        finalsequenzaoutput="analysis/clonality/{run}/{run}.bin50.final.seqz.txt.gz"
    benchmark:
        "benchmarks/clonality/{run}/{run}_clonality_addheader.txt"
    shell:
        #"zcat {input.headerfile} {input.binned_file_in} |bgzip > {output.finalsequenzaoutput}"
         "cat {input.headerfile} {input.binned_file_in} > {output.finalsequenzaoutput}"

rule sequenza_fileprep:
    input:
        bin50="analysis/clonality/{run}/{run}.bin50.final.seqz.txt.gz"
    params:
        out_dir="analysis/clonality",
        sample_name=lambda wildcards:wildcards.run,
        sequenza_bin_path="~/miniconda3/envs/sequenza/bin/",
    output:
        pyclone_tsv="analysis/clonality/{run}/{run}_pyclone.tsv",
        genome_view_plot="analysis/clonality/{run}/{run}_genome_view.pdf",
    benchmark:
        "benchmarks/clonality/{run}_{run}_postRscript.txt"
    shell:
        #"""{params.sequenza_bin_path}Rscript cidc_wes/modules/scripts/sequenza.R  {input}  {params.out_dir}/{wildcards.run}  {params.sample_name}"""
        """{params.sequenza_bin_path}Rscript /mnt/ssd/mda-r1-pt1_report/cidc_wes/modules/scripts/sequenza.R  {input}  {params.out_dir}/{wildcards.run}  {params.sample_name}"""

rule pyclone_build_mutation_file:
    input:
        "analysis/clonality/{run}/{run}_pyclone.tsv"
    params:
        pyclone_bin_path="~/miniconda3/envs/pyclone/bin/",
    output:
        "analysis/clonality/{run}/{run}_pyclone.yaml"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone_build_mutation_file.txt"
    shell:
        "{params.pyclone_bin_path}PyClone build_mutations_file --in_file {input} --out_file {output}"

rule pyclone_generate_config:
    input:
        "analysis/clonality/{run}/{run}_pyclone.yaml"
    output:
        "analysis/clonality/{run}/pyclone.config.yaml"
    params:
        outdir = lambda wildcards: "analysis/clonality/%s" % wildcards.run,
        #HARD-coding
        template = "/mnt/ssd/mda-r1-pt1_report/cidc_wes/static/clonality/pyclone.config.yaml"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone.generate.config.txt"
    shell:
        "/mnt/ssd/mda-r1-pt1_report/cidc_wes/modules/scripts/clonality_pycloneConfig.py -t {params.template} -m {input} -o {params.outdir} > {output}"

rule pyclone_run_analysis:
    input:
        "analysis/clonality/{run}/pyclone.config.yaml"
    params:
        pyclone_bin_path="~/miniconda3/envs/pyclone/bin/",
    output:
        #LEN: I don't know what to key in on as output
        "analysis/clonality/{run}/trace/alpha.tsv.bz2"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone.run_analysis.txt"
    shell:
        "{params.pyclone_bin_path}PyClone run_analysis --config_file {input}"
