#module: somatic tumor calling using TNscope

def somatic_calling_tumor_TNscope_inputFn(wildcards):
    run = config['runs'][wildcards.run]
    tumor = run[1]
    ls = ["analysis/align/%s/%s_recalibrated.bam" % (tumor, tumor)]
    return ls

def somatic_calling_tumor_TNscope_filter_inputFn(wildcards):
    run = wildcards.run
    ls = ["analysis/somatic/%s/%s_preprocess.vcf.gz" % (run, run)]
    return ls


if config.get('tumor_only'):
    rule somatic_calling_tumor_TNscope:
        input:
            realignedbam = somatic_calling_tumor_TNscope_inputFn
        output:
            tnscopevcf="analysis/somatic/{run}/{run}_preprocess.vcf.gz"
        params:
            index=config['genome_fasta'],
            sentieon_path=config['sentieon_path'],
	        pon_haplotyper= config['pons_cidc'],
            dbsnp= config['dbsnp'],
            tumor = lambda wildcards: config['runs'][wildcards.run][1],
            trim_soft_clip = "--trim_soft_clip" if config.get("trim_soft_clip", False) else "",
        threads: 18 #_somatic_threads
        priority: 50
        group: "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}.somatic_calling_TNscope.txt"
        shell:
            """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input} --algo TNscope --tumor_sample {params.tumor} --pon {params.pon_haplotyper}  --dbsnp {params.dbsnp} {params.trim_soft_clip} {output.tnscopevcf}"""
    rule somatic_calling_tumor_filter:
        input:
            somatic_calling_tumor_TNscope_filter_inputFn
        output:
            vcf= "analysis/somatic/{run}/{run}_filtered.vcf",
            gz_vcf= "analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
        params:
            pon_gnomad = config['gnomad'],
            outdir = "analysis/somatic/{run}"
        group: "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}.somatic_calling_filter_TNscope.txt"
        shell:
            "bcftools isec -C -p {params.outdir} {input}  {params.pon_gnomad}"
            "&& cp  {params.outdir}/0000.vcf {output.vcf}"
            "&& bgzip -c {output.vcf} > {output.gz_vcf}"

    rule somatic_creating_tbi:
        input:
            "analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
        output:
            "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz.tbi"
        group:
            "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}_{caller}.somatic_creating_tbi.txt"
        shell:
            "tabix -p vcf {input}"



else: #Except tumor-normal pairing for samples
    rule somatic_calling_tumor_TNscope:
        input:
            corealignedbam="analysis/corealignments/{run}/{run}_tn_corealigned.bam"
        output:
            tnscopevcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz",
	    tnscopevcf_tbi="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz.tbi"
        params:
            index=config['genome_fasta'],
            sentieon_path=config['sentieon_path'],
            dbsnp= config['dbsnp'],
            tumor = lambda wildcards: config['runs'][wildcards.run][1],
            normal = lambda wildcards: config['runs'][wildcards.run][0],
            trim_soft_clip = "--trim_soft_clip" if config.get("trim_soft_clip", False) else "",
        threads: 18 #_somatic_threads
        priority: 50
        group: "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}.somatic_calling_TNscope.txt"
        shell:
            """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input.corealignedbam} --algo TNscope --tumor_sample {params.tumor} --normal_sample {params.normal} --dbsnp {params.dbsnp} {params.trim_soft_clip} {output.tnscopevcf}"""
