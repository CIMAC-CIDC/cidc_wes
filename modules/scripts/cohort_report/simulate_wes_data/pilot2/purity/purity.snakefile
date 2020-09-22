def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("json/%s.purity.json" % run)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule maf_json:
    input:
        "raw_data/{run}.optimalpurityvalue.txt"
    output:
        "json/{run}.purity.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./purity_json_writer.py -r {params.run} -f {input} -o {output}"
