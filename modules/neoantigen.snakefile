#module: neoantigen prediction module using pvacseq pipeline (from pvactools)

_neoantigen_threads=8 #need to tune this

def neoantigen_runsHelper(wildcards, iindex):
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


def neoantigen_getNormal(wildcards):
    return neoantigen_runsHelper(wildcards, 0)

def neoantigen_getTumor(wildcards):
    return neoantigen_runsHelper(wildcards, 1)

def neoantigen_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        tumor = config['runs'][run][1]
        ls.append("analysis/neoantigen/%s/MHC_Class_I/%s.filtered.condensed.ranked.tsv" % (run,tumor))
    return ls

def getPvacseqOut(wildcards):
    """returns tuple (run, tumor sample name)"""
    run = wildcards.run
    tumor = config['runs'][run][1]
    return(run,tumor)

rule neoantigen_all:
    input:
        neoantigen_targets

rule neoantigen_vep_annotate:
    input:
        "analysis/somaticVariants/{run}/{run}_tnsnv.filter.vcf"
    output:
        "analysis/somaticVariants/{run}/{run}_tnsnv.filter.neoantigen.vep.vcf"
    params:
        index=config['genome_fasta'],
        vep_data=config['vep_data'],
        vep_plugins=config['vep_plugins'],

        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
    benchmark:
        "benchmarks/neoantigen/{run}/{run}.neoantigen_vep_annotate.txt"
    shell:
        """vep --input_file {input} --output_file {output} --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta {params.index} --offline --cache --dir_cache {params.vep_data} --plugin Downstream --plugin Wildtype --dir_plugins {params.vep_plugins} --pick"""


rule neoantigen_pvacseq:
    input:
        "analysis/somaticVariants/{run}/{run}_tnsnv.filter.neoantigen.vep.vcf"
    output:
        "analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.condensed.ranked.tsv"
        #NOTE: the wildcard tumor is used b/c we can't actually pull in the
        #tumor name 

        #NOTE: typically when HLA class I and HLA class II are both called
        #the results are in combined/{inputname}.condensed.tsv
        #Since we're only generating class I for now, no combined is generated
    params:
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
        HLA = "HLA-A*02:01,HLA-B*35:01,DRB1*11:01", #Hard coding antigents
        callers=" MHCflurry NetMHC", #hardcoding list of callers to try
        output_dir = lambda wildcards: "analysis/neoantigen/%s/" % wildcards.run,
        epitope_lengths="8,9,10", #hard-coding these--these are standard
    threads: _neoantigen_threads
    benchmark:
        "benchmarks/neoantigen/{run}/{run}.neoantigen_pvacseq.txt"
    shell:
        """pvacseq run {input} {params.tumor} {params.HLA} {params.callers} {params.output_dir} -e {params.epitope_lengths} -t {threads} --normal-sample-name {params.normal}"""
  
