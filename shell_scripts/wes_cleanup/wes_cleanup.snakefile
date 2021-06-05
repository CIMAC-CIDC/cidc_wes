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

rule cleanupBucket:
    output:
        "analysis/cleanup/{sample}.txt"
    params:
        bucket_path = lambda wildcards: getBucketPath(wildcards)
    shell:
        """./wes_cleanup.py -b {params.bucket_path} &&
        touch {output}"""
    
