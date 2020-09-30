def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("json/%s.hla.json" % run)
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
    tmp['tumor_opti'] = "analysis/optitype/%s/%s_result.tsv" % (tmr, tmr)
    tmp['tumor_xhla'] = "analysis/xhla/%s/report-%s-hla.json" % (tmr, tmr)

    tmp['normal_opti'] = "analysis/optitype/%s/%s_result.tsv" % (nrm, nrm)
    tmp['normal_xhla'] = "analysis/xhla/%s/report-%s-hla.json" % (nrm, nrm)

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
        "json/{run}.hla.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./hla_json_writer.py -r {params.run} -t {input.tumor_opti} -s {input.tumor_xhla} -n {input.normal_opti} -m {input.normal_xhla} -o {output}"
