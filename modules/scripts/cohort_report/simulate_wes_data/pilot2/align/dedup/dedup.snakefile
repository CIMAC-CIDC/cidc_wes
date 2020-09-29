def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        tmr = config[run]['tumor']
        nrm = config[run]['normal']
        ls.append("mapping/%s_dedup.txt" % tmr)
        ls.append("mapping/%s_dedup.txt" % nrm)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule map_stats:
    input:
        "raw_data/{sample}.sorted.dedup.bam"
    output:
        "mapping/{sample}_dedup.txt"
    threads: 8
    shell:
        "sambamba flagstat -t {threads} {input} > {output}"
