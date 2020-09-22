def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("json/%s.filtered_maf.json" % run)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

#NOT needed for run-level
# def getFiles(wildcards):
#     """Will return a list of the run's tumor and normal files"""
#     tmp = {}
#     run = config[wildcards.run]
#     tmp['tumor'] = "raw_data/%s_coverage_metrics.txt.sample_summary" % run['tumor']
#     tmp['normal'] = "raw_data/%s_coverage_metrics.txt.sample_summary" % run['normal']
#     #print(tmp)
#     return tmp

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule maf_json:
    input:
        "raw_data/{run}_tnscope.filter.maf"
    output:
        "json/{run}.filtered_maf.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./filtered_maf_json_writer.py -r {params.run} -m {input} -o {output}"
