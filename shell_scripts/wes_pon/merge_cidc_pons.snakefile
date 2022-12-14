# Len Taing 2021 (TGBTG)
# Generate a panel of normal according to
# ref: https://support.sentieon.com/manual/TNscope_usage/tnscope/

configfile: "merge_cidc_pons.config.yaml"
def targets(wildcards):
    ls = []
    ls.append("analysis/pon/cidc.pon.vcf.gz")
    ls.append("analysis/pon/cidc.pon.vcf.gz.tbi")
    return ls

rule all:
    input:
        targets

rule mergeVCFs:
    input:
        expand("{pon_file}", pon_file=config['pon_files'])
    output:
        "analysis/pon/{center}.pon.vcf"
    benchmark:
        "benchmarks/mergeVCFs.{center}.txt"
    shell:
        "bcftools merge -m all -f PASS,. --force-samples {input} | "
        "bcftools plugin fill-AN-AC | "
        "bcftools filter -i 'SUM(AC)>1' > {output}"

rule gzipAndTabix:
    input:
        "analysis/pon/{center}.pon.vcf"
    output:
        gz="analysis/pon/{center}.pon.vcf.gz",
        tbi="analysis/pon/{center}.pon.vcf.gz.tbi",
    benchmark:
        "benchmarks/gzipAndTabix.{center}.txt"
    shell:
        #"bgzip -c {input} > {output.gz} && tabix -p vcf {output.gz}"
        "bgzip {input} && tabix -p vcf {output.gz}"
