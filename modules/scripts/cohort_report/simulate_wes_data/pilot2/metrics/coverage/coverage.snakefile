def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("json/%s.coverage.json" % run)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

def getFiles(wildcards):
    """Will return a list of the run's tumor and normal files"""
    tmp = {}
    run = config[wildcards.run]
    tmp['tumor'] = "raw_data/%s_coverage_metrics.txt.sample_summary" % run['tumor']
    tmp['normal'] = "raw_data/%s_coverage_metrics.txt.sample_summary" % run['normal']
    #print(tmp)
    return tmp

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule gc_content_json:
    input:
        unpack(getFiles)
    output:
        "json/{run}.coverage.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./coverage_json_writer.py -r {params.run} -t {input.tumor} -n {input.normal} -o {output}"
