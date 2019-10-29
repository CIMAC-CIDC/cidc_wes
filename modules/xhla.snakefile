# module:Precision HLA typing from next-generation sequencing data by xHLA
_xhla_threads=16

def xhla_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/xhla/%s/report-%s-hla.json" % (sample,sample))
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
        chr6sortbamfile = "analysis/xhla/{sample}/report-{sample}-hla.json"
    threads:_xhla_threads
    group: "xhla"
    params:
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "%sanalysis/xhla/%s/" % (config['remote_path'], wildcards.sample),
        path="source activate %s" % config['xhla_root'],
    #singularity: "docker://humanlongevity/hla"
    conda: "../envs/xhla_env.yml"
    benchmark:
        "benchmarks/xhla/{sample}/{sample}.xhla.txt"
    shell:
        """{params.path}; run.py --sample_id {params.name} --input_bam_path {input.in_sortbamfile} --output_path  {params.output_dir}"""
