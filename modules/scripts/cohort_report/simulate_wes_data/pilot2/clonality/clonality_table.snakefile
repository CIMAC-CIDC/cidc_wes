def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/clonality/%s/%s_table.tsv" % (run,run))
        #tmr = config[run]['tumor']
        #nrm = config[run]['normal']
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (tmr,tmr))
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.bam" % (nrm,nrm))
        
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule pyclone_build_table:
    input:
        conf="analysis/clonality/{run}/pyclone.config.yaml",
        other="analysis/clonality/{run}/trace/alpha.tsv.bz2",
        #other="raw_data/{run}/alpha.tsv.bz2",
    params:
        pyclone_bin_path="~/miniconda3/envs/pyclone/bin/",
    output:
        "analysis/clonality/{run}/{run}_table.tsv"
    benchmark:
        "benchmarks/clonality/{run}/{run}.pyclone.build_table.txt"
    shell:
        "{params.pyclone_bin_path}PyClone build_table --config_file {input.conf} --out_file {output} --table_type cluster --max_clusters 100"
