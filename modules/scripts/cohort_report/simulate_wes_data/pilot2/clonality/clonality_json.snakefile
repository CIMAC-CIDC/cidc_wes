def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("json/%s.clonality.json" % run)
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
        "analysis/clonality/{run}/{run}_table.tsv"
    output:
        "json/{run}.clonality.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./clonality_json_writer.py -r {params.run} -f {input} -o {output}"
