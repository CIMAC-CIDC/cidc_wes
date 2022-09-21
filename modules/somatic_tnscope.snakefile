#module: somatic tumor calling using TNscope

def somatic_calling_tumor_TNscope_inputFn(wildcards):
    run = config['runs'][wildcards.run]
    tumor = run[1]
    ls = ["analysis/align/%s/%s_recalibrated.bam" % (tumor, tumor)]
    return ls

#DEPRECATED
#def somatic_calling_tumor_TNscope_filter_inputFn(wildcards):
#    run = wildcards.run
#    ls = ["analysis/somatic/%s/%s_preprocess.vcf.gz" % (run, run)]
#    return ls


if config.get('tumor_only'):
    rule somatic_calling_tumor_TNscope:
        input:
            realignedbam = somatic_calling_tumor_TNscope_inputFn
        output:
            tnscopevcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz",
        params:
            index=config['genome_fasta'],
            sentieon_path=config['sentieon_path'],
	    pon= config['pons_cidc'],
            dbsnp= config['dbsnp'],
            tumor = lambda wildcards: config['runs'][wildcards.run][1],
            trim_soft_clip = "--trim_soft_clip" if config.get("trim_soft_clip", False) else "",
        threads: 18 #_somatic_threads
        priority: 50
        group: "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}.somatic_calling_TNscope.txt"
        shell:
            """{params.sentieon_path}/sentieon driver -r {params.index} -t {threads}  -i {input} --algo TNscope --tumor_sample {params.tumor} --pon {params.pon}  --dbsnp {params.dbsnp} {params.trim_soft_clip} {output.tnscopevcf}"""
    rule somatic_tumor_only_filters:
        """filters out variants in panel_of_normal and in dbSNP ('DB;')"""
        input:
            "analysis/somatic/{run}/{run}_tnscope.output.vcf.gz"
        output:
            "analysis/somatic/{run}/{run}_preprocess.vcf"
        group: "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}.somatic_tumor_only_filters"
        shell:
            "zgrep -v 'panel_of_normals\|DB;' {input} > {output}"

    rule somatic_bgzip_tbi:
        "bgzip and index {run}_preprocess.vcf"
        input:
            "analysis/somatic/{run}/{run}_preprocess.vcf"
        output:
            vcf="analysis/somatic/{run}/{run}_preprocess.vcf.gz",
            tbi="analysis/somatic/{run}/{run}_preprocess.vcf.gz.tbi",
        group:
            "somatic"
        benchmark:
            "benchmarks/somatic/{run}/{run}.somatic_bgzip_tbi.txt"
        shell:
            "bgzip -c {input} > {output.vcf} && tabix -p vcf {output.vcf}"

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
