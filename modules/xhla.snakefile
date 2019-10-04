# module:Precision HLA typing from next-generation sequencing data by xHLA
_xhla_threads=16

def xhla_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/xhla/hla-%s/%s.json" % (sample,sample))
        #ls.append("analysis/optitype/hla-%s/%s.tsv" % (sample,sample))
        #ls.append("analysis/optitype/hla-%s/%s.hla" % (sample,sample))
        #ls.append("analysis/optitype/hla-%s/%s.fq.gz" % (sample,sample))
        #ls.append("analysis/optitype/%s/%s.sorted.chr6.end2.fastq" % (sample,sample))
    return ls

rule xhla_all:
    input:
        xhla_targets

rule xhla:
    """calculate hlatyping by xhla"""
    input:
        in_sortbamfile = "analysis/align/{sample}/{sample}.sorted.bam"
    output:
        chr6sortbamfile = "analysis/xhla/hla-{sample}/{sample}.json"
    threads:_xhla_threads
    group: "xhla"
    params:
        name=lambda wildcards: wildcards.sample,
        #output_dir=lambda wildcards: "%sanalysis/xhla/%s/" % (config['remote_path'], wildcards.sample),
        xhla_path=config['xhla_path']
    singularity: "docker://humanlongevity/hla"
    benchmark:
        "benchmarks/xhla/{sample}/{sample}.xhla.txt"
    shell:
        """{params.xhla_path}/typer.sh {input.in_sortbamfile} {params.name}"""
