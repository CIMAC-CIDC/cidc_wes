#MODULE: wes report2 module 
from yaml import dump as yaml_dump

def report2_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #META
    ls.append("analysis/report2/wes_meta/wes_run_version.tsv")
    ls.append("analysis/report2/wes_meta/wes_software_versions.tsv")
    ls.append("analysis/report2/wes_meta/wes_meta.yaml")
    #ALIGN
    ls.append("analysis/report2/alignment/mapping_stats.tsv")    
    return ls

rule report2_all:
    input:
        "analysis/report2/report.html"

###############################################################################
#META information
rule report2_meta_version:
    """Gather's all the information required for the meta information and
    outputs a csv file
    """
    input:
        config="config.yaml",
        wes_versions="cidc_wes/static/wes_versions.yaml"
    output:
         "analysis/report2/wes_meta/wes_run_version.tsv"
    message:
        "REPORT: creating WES version table"
    group: "report2"
    shell:
        """cidc_wes/modules/scripts/report_meta.py -c {input.config} -v {input.wes_versions} -o {output}"""

rule report2_meta_software:
    """Gather's all the information required for the meta information and
    outputs a csv file
    """
    input:
        config="config.yaml",
        wes_versions="cidc_wes/static/wes_versions.yaml"
    output:
         "analysis/report2/wes_meta/wes_software_versions.tsv"
    message:
        "REPORT: creating WES software table"
    group: "report2"
    shell:
        """cidc_wes/modules/scripts/report_software.py -c {input.config} -v {input.wes_versions} -o {output}"""

rule report2_meta_order:
    """Dictate the ordering of these parts"""
    input:
        "analysis/report2/wes_meta/wes_software_versions.tsv",
        "analysis/report2/wes_meta/wes_run_version.tsv"
    params:
        #order= yaml_dump({'order': ['wes_run_version.tsv','wes_software_versions.tsv']})
        order= yaml_dump({'order': ['wes_software_versions.tsv', 'wes_run_version.tsv']})
    output:
        "analysis/report2/wes_meta/wes_meta.yaml"
    group: "report2"
    shell:
        """echo '{params.order}' > {output}"""
###############################################################################

rule report2_alignment_table:
    """Generate the mapping stats table for the report"""
    input:
        "analysis/align/mapping.csv"
    output:
        tsv="analysis/report2/alignment/mapping_stats.tsv",
        cap="analysis/report2/alignment/mapping_stats_caption.txt",
    params:
        caption="This table shows the Total number reads in each sample and how many of those reads were mapped."
    message:
        "REPORT: creating mapping stats for alignment section"
    group: "report2"
    shell:
        """echo "{params.caption}" > {output.cap} && cidc_wes/modules/scripts/report_align_mappingStats.py -f {input} -o {output.tsv}"""

###############################################################################
rule report2_slurp:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report2_targets
    params:
        report_path = "analysis/report2",
        sections_list=",".join(['wes_meta','alignment'])
    output:
        "analysis/report2/report.html"
    message:
        "REPORT: Generating WES report"
    group: "report2"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -o {output} && cp -r cidc_wes/report2/static {params.report_path}"""

