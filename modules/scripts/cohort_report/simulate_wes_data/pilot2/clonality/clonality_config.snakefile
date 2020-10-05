def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/pyclone.config.yaml" % run)
        #tmr = config[run]['tumor']
        #nrm = config[run]['normal']
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (tmr,tmr))
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (nrm,nrm))
        
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

configfile: "config.cohort_json.big.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule pyclone_build_mutation_file:
    input:
        "tsv/{run}_pyclone.tsv"
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
        template = "/mnt/ssd/mda-r1-pt1_report/cidc_wes/static/clonality/pyclone.config.yaml"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone.generate.config.txt"
    shell:
        "/mnt/ssd/mda-r1-pt1_report/cidc_wes/modules/scripts/clonality_pycloneConfig.py -t {params.template} -m {input} -o {params.outdir} > {output}"

