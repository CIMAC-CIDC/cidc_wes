def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("json/%s.is_metrics.json" % run)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

def getISFiles(wildcards):
    """Will return a list of the run's tumor and normal files"""
    tmp = {}
    run = config[wildcards.run]
    tmp['tumor'] = "raw_data/%s_is_metrics.txt" % run['tumor']
    tmp['normal'] = "raw_data/%s_is_metrics.txt" % run['normal']
    #print(tmp)
    return tmp

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule gc_content_json:
    input:
        unpack(getISFiles)
    output:
        "json/{run}.is_metrics.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./insert_size_json_writer.py -r {params.run} -t {input.tumor} -n {input.normal} -o {output}"
