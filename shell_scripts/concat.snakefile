configfile: "10104_concat.config.yaml"

from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

def targets(wildcards):
    ls = []
    for s in config['samples']:
        ls.append(GS.remote("%s/%s.fastq.gz" % (config["destination"],s)))
    return ls

rule all:
    input:
        targets

def concatInput(wildcards):
    s = wildcards.sample
    dict={}
    dict['L1']=GS.remote(config['samples'][s][0])
    dict['L2']=GS.remote(config['samples'][s][1])
    return dict

rule concat:
    input:
        unpack(concatInput)
    output:
        GS.remote("%s/{sample,[^/]+}.fastq.gz" % config["destination"])
    shell:
        "cat {input.L1} {input.L2} > {output}"

