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
    tmr = run['tumor']
    nrm = run['normal']
    tmp['tumor_map'] = "mapping/raw_data/%s_mapping.txt" % tmr
    tmp['tumor_dedup'] = "dedup/mapping/%s_dedup.txt" % tmr

    tmp['normal_map'] = "mapping/raw_data/%s_mapping.txt" % nrm
    tmp['normal_dedup'] = "dedup/mapping/%s_dedup.txt" % nrm

    #print(tmp)
    return tmp

configfile: "config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule mapping_json:
    input:
        unpack(getFiles)
    output:
        "json/{run}.coverage.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./mapping_json_writer.py -r {params.run} -t {input.tumor_map} -s {input.tumor_dedup} -n {input.normal_map} -m {input.normal_dedup} -o {output}"
