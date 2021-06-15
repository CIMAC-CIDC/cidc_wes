configfile: "config.wes_cleanup.yaml"

def targets(wildcards):
    ls = []
    for sample in config['samples']:
        ls.append("analysis/cleanup/%s.txt" % sample)
    return ls

rule all:
    input:
        targets

def getBucketPath(wildcards):
    sample = wildcards.sample
    bp = config['samples'][sample][0]
    #print(bp)
    return bp

def isTumorOnly(wildcards):
    ret = ""
    tumor_only = config.get('tumor_only', [])
    if tumor_only and wildcards.sample in tumor_only:
        ret = " -t "
    return ret

rule cleanupBucket:
    output:
        "analysis/cleanup/{sample}.txt"
    params:
        bucket_path = lambda wildcards: getBucketPath(wildcards),
        tumor_only = lambda wildcards: isTumorOnly(wildcards)
    shell:
        """./wes_cleanup.py -b {params.bucket_path} {params.tumor_only} &&
        touch {output}"""
    
