rule somatic_calling_TNhaplotyper2:
    input:
        tumorbam=somatic_getNTumor_recal,
	tumorbai=somatic_getNTumor_recal_bai,
        normalbam=somatic_getNormal_recal,
        normalbai=somatic_getNormal_recal_bai,
    output:
        tnhaplotyper2vcf="analysis/somatic/{run}/{run}_tnhaplotyper2.output.vcf.gz"
        #LEN/AASHNA: IF tnhaplotyper2 generates tbi, we should add it as output
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        tnhaplotyper_pon= config['pons_haplotyper'],
	#JUST sample names - can also use the helper fns, e.g.
	#normal = lambda wildcards: getNormal_sample(wildcards)
	normal = lambda wildcards: config['runs'][wildcards.run][0],
	tumor = lambda wildcards: config['runs'][wildcards.run][1],
        trim_soft_clip = "--trim_soft_clip" if config.get("trim_soft_clip", False) else "",
    threads:_somatic_threads
    benchmark:
        "benchmarks/somatic/{run}/{run}.somatic_calling_TNhaplotyper2.txt"
    shell:
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index}  -i {input.normalbam}  -i  {input.tumorbam}  --algo TNhaplotyper2 --pon {params.tnhaplotyper_pon}  --tumor_sample {params.tumor} --normal_sample {params.normal} {params.trim_soft_clip} {output.tnhaplotyper2vcf}"""
