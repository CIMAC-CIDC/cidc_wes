rule somatic_calling_TNscope:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        tnscopevcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        #JUST sample names - can also use the helper fns, e.g.
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
        trim_soft_clip = "--trim_soft_clip" if config.get("trim_soft_clip", False) else "",
    threads: 18 #_somatic_threads
    priority: 50
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}.somatic_calling_TNscope.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNscope --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {params.trim_soft_clip} {output.tnscopevcf}"""
