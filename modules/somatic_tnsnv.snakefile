rule somatic_calling_TNsnv:
    input:
        corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
    output:
        statscall="analysis/somatic/{run}/{run}_call.output.stats",
        tnsnvvcf="analysis/somatic/{run}/{run}_tnsnv.output.vcf.gz",
        tnsnvvcftbi="analysis/somatic/{run}/{run}_tnsnv.output.vcf.gz.tbi",
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        dbsnp= config['dbsnp'],
        #JUST sample names - can also use the helper fns, e.g.
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
    threads:96
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}.somatic_calling_TNsnv.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads} -i {input.corealignedbam} --algo TNsnv --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} --call_stats_out {output.statscall} {output.tnsnvvcf}"""
