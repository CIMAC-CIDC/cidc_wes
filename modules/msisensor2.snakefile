# author: Len Taing (TGBTG) 
# year: 2020
# module: microsatellite instablility analysis using msisensor2

#Deprecated
# def msisensor2_getNormalInput(wildcards):
#     run = config['runs'][wildcards.run]
#     normal = run[0]
#     return "analysis/align/%s/%s.sorted.dedup.bam" % (normal, normal)

# def msisensor2_getTumorInput(wildcards):
#     run = config['runs'][wildcards.run]
#     tumor = run[1]
#     return "analysis/align/%s/%s.sorted.dedup.bam" % (tumor, tumor)

def msisensor2_output_files(wildcards):
    """returns a list of filepaths generated by this module to store 
    in the CIDC for a given sample 
    """
    ls = []
    run = wildcards.run
    ls.append("analysis/msisensor2/%s/%s_msisensor" % (run,run))
    return ls

def msisensor2_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        ls.append("analysis/msisensor2/%s/%s_msisensor.txt" % (run,run))
        #output file map
        ls.append("analysis/msisensor2/%s/%s.msisensor2.output.yaml" % (run,run))
    return ls

rule msisensor2_all:
    input:
        msisensor2_targets
    benchmark: "benchmarks/msisensor2/msisensor2_all.txt"

def msisensor2_msisensor2InputFn(wildcards):
    run = config['runs'][wildcards.run]
    normal = run[0]
    tumor = run[1]

    tmp = {}
    if not config.get('tumor_only'): #Only run when we have normals
        tmp['normal'] = "analysis/align/%s/%s.sorted.dedup.bam" % (normal, normal)
    tmp['tumor'] = "analysis/align/%s/%s.sorted.dedup.bam" % (tumor, tumor)
    return tmp

rule msisensor2:
    """calculate microsatellite instability with msisensor2"""
    input:
        unpack(msisensor2_msisensor2InputFn)
    output:
        "analysis/msisensor2/{run}/{run}_msisensor",
        "analysis/msisensor2/{run}/{run}_msisensor_dis",
        "analysis/msisensor2/{run}/{run}_msisensor_somatic",
        #"analysis/msisensor2/{run}/{run}_msisensor_germline",
    group: "msisensor2"
    params:
        #Change the files inputted depending on whether normal is available 
        in_files = lambda wildcards,input: "-t %s" % input.tumor if config.get('tumor_only', False) else "-t %s -n %s" % (input.tumor, input.normal),
        outpath = "analysis/msisensor2/{run}/{run}_msisensor",
        ref_path= config["msisensor2"],
    #conda: "../envs/msisensor2_env.yml"
    log: "analysis/logs/msisensor2/{run}/{run}_msisensor2.log.txt"
    benchmark:
        "benchmarks/msisensor2/{run}/{run}.msisensor2.txt"
    shell:
        """msisensor2 msi -M {params.ref_path} {params.in_files} -o {params.outpath} > {log}"""

rule msisensor2_make_file_map:
    input:
        msisensor2_output_files
    output:
        "analysis/msisensor2/{run}/{run}.msisensor2.output.yaml"
    benchmark: "benchmarks/msisensor2/{run}/{run}.msisensor2_make_file_map.txt"
    group: "msisensor2"
    params:
        run = lambda wildcards: wildcards.run,
        kkeys = " -k ".join(['msisensor2_results']),
        files = lambda wildcards, input: " -f ".join(input),
    shell:
        "cidc_wes/modules/scripts/yaml_writer.py -t runs -n {params.run} -k {params.kkeys} -f {params.files} > {output}"

rule msisensor2_copy:
    """Rename analysis/msisensor2/{run}/{run}_msisensor as 
    analysis/msisensor2/{run}/{run}_msisensor.txt"""
    input:
        "analysis/msisensor2/{run}/{run}_msisensor",
    output:
        "analysis/msisensor2/{run}/{run}_msisensor.txt",
    group: "msisensor2"
    shell:
        "cp {input} {output}"
