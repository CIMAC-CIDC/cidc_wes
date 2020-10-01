def targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("cohort/%s.json" % run)
    return ls

def getRuns(config):
    """Will return a list of samples"""
    return list(config.keys())

def getFiles(wildcards):
    """Will return a list of the run's json files"""
    run = wildcards.run #config[wildcards.run]
    #tmr = run['tumor']
    #nrm = run['normal']
    tmp = {}
    tmp['mapping'] = "align/json/%s.mapping.json" % run
    tmp['coverage'] = "metrics/coverage/json/%s.coverage.json" % run
    tmp['gc'] = "metrics/gc_metrics/json/%s.gc_metrics.json" % run
    tmp['is'] = "metrics/is_metrics/json/%s.is_metrics.json" % run
    tmp['mean_quality']="metrics/mean_quality/json/%s.mean_quality.json" % run
    tmp['hla']="hla/json/%s.hla.json" % run
    tmp['purity']="purity/json/%s.purity.json" % run
    tmp['somatic']="somatic/maf/json/%s.filtered_maf.json" % run
    #FAKE clonality
    tmp['clonality']="clonality/clonality_stub.json"

    #print(tmp)
    return tmp

configfile: "align/config.cohort_json.yaml"
config['runs'] = getRuns(config)

rule all:
    input:
        targets

rule stitch_json:
    input:
        unpack(getFiles)
    output:
        "cohort/{run}.json"
    params:
        run = lambda wildcards: wildcards.run
    shell:
        "./stitcher_json_writer.py -r {params.run} -m {input.mapping} -c {input.coverage} -g {input.gc} -i {input.is} -q {input.mean_quality} -j {input.hla} -p {input.purity} -s {input.somatic} -t {input.clonality} -o {output}"
