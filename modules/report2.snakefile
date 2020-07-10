#MODULE: wes report2 module 
from yaml import dump as yaml_dump

def report2_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #META
    ls.append("analysis/report2/wes_meta/wes_run_version.tsv")
    ls.append("analysis/report2/wes_meta/wes_software_versions.tsv")
    ls.append("analysis/report2/wes_meta/wes_meta.yaml")
    #Data Quality
    ls.append("analysis/report2/data_quality/mapping_stats.tsv")
    ls.append("analysis/report2/data_quality/coverage_statistics.tsv")
    ls.append("analysis/report2/data_quality/data_quality.yaml")
    ls.append("analysis/report2/somatic_variants/summary_table.csv")
    ls.append("analysis/report2/somatic_variants/functional_annotation.csv")
    ls.append("analysis/report2/somatic_variants/SNV_statistics.csv")
    ls.append("analysis/report2/somatic_variants/somatic_variants.yaml")
    ls.append("analysis/report2/somatic_variants/tumor_mutational_burden.tsv")

    ls.append("analysis/report2/neoantigens/HLA_results.tsv")
    ls.append("analysis/report2/neoantigens/neoantigen_list.tsv")

    ls.append("analysis/report2/copy_number/tumor_purity.tsv")
    ls.append("analysis/report2/copy_number/tumor_clonality.tsv")
    for run in config['runs']:
        ls.append("analysis/report2/data_quality/%s-qc_plots.tsv" % run)
        ls.append("analysis/report2/somatic_variants/%s_lego_plot.png" % run)
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

rule report2_data_quality_table:
    """Generate the mapping stats table for the report"""
    input:
        "analysis/align/mapping.csv"
    output:
        tsv="analysis/report2/data_quality/mapping_stats.tsv",
        cap="analysis/report2/data_quality/mapping_stats_caption.txt",
    params:
        caption="This table shows the Total number reads in each sample and how many of those reads were mapped."
    message:
        "REPORT: creating mapping stats for data_quality section"
    group: "report2"
    shell:
        """echo "{params.caption}" > {output.cap} && cidc_wes/modules/scripts/report_dataQual_mappingStats.py -f {input} -o {output.tsv}"""

#------------------------------------------------------------------------------
def report2_data_quality_plotsInputFn(wildcards):
    """Given a run, returns a list of the various sub plots to generate
    """
    ret = []
    run = wildcards.run
    for sample in config['runs'][wildcards.run]:
        ret.append("analysis/report2/data_quality/%s-qc_plots/%s_gcBias.png" % (run, sample))
        ret.append("analysis/report2/data_quality/%s-qc_plots/%s_qualityScore.png" % (run, sample))
        ret.append("analysis/report2/data_quality/%s-qc_plots/%s_qualityByCycle.png" % (run, sample))
        ret.append("analysis/report2/data_quality/%s-qc_plots/%s_insertSize.png" % (run, sample))
    return ret

rule report2_data_quality_gcPlot:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report2/data_quality/{run}-qc_plots/{sample}_gcBias.png"
    params:
        page = 1
    message:
        "REPORT: generating data_quality gc bias plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report2_data_quality_qualityScore:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report2/data_quality/{run}-qc_plots/{sample}_qualityScore.png"
    params:
        page = 2
    message:
        "REPORT: generating data_quality quality score plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report2_data_quality_qualityByCycle:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report2/data_quality/{run}-qc_plots/{sample}_qualityByCycle.png"
    params:
        page = 3
    message:
        "REPORT: generating data_quality quality by cycle plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report2_data_quality_insertSize:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report2/data_quality/{run}-qc_plots/{sample}_insertSize.png"
    params:
        page = 4
    message:
        "REPORT: generating data_quality insert size plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report2_data_quality_plots_table:
    """Generate the fastqc plots table for the report
    The trick here is that the plot table shows BOTH samples, tumor and 
    normal, for the run--so we need 1. an inpput fn to look up the samples
    and 2. params to generate the correct plots
    """
    input:
        report2_data_quality_plotsInputFn
    output:
        tsv="analysis/report2/data_quality/{run}-qc_plots.tsv",
        sub_cap = "analysis/report2/data_quality/{run}-qc_plots_subcaption.txt",
        #cap=...
    params:
        normal= lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
        image_paths = lambda wildcards: "analysis/report2/data_quality/%s-qc_plots/" % (wildcards.run),
        sub_caption = "NOTE: (T) denotes tumor sample; (N) denotes normal sample",
    message:
        "REPORT: creating QC plots for data_quality section"
    group: "report2"
    shell:
        """echo "{params.sub_caption}" > {output.sub_cap} && cidc_wes/modules/scripts/report_dataQual_plot_table.py -n {params.normal} -t {params.tumor} -p {params.image_paths} -o {output}"""

rule report2_data_quality_coverage:
    """Generate the coverage data table"""
    input:
        "analysis/metrics/all_sample_summaries.txt"
    output:
        tsv="analysis/report2/data_quality/coverage_statistics.tsv",
        cap="analysis/report2/data_quality/coverage_statistics_caption.txt"
    params:
        caption="""The following table describes the read depth coverage statistics. With the exception of the "Total Reads" column, which represents the total number of reads in each sample, all numbers represent reads in targeted regions."""
    shell:
        """echo "{params.caption}" > {output.cap} && cidc_wes/modules/scripts/report_dataQual_coverage.py -f {input} -o {output.tsv}"""
    
def report2_data_quality_orderInputFn(wildcards):
    """Returns the first run in config 
    NOTE: we're running WES with one run at a time
    """
    ret = ["analysis/report2/data_quality/mapping_stats.tsv"]
    run = list(config['runs'].keys())[0]
    ret.append("analysis/report2/data_quality/%s-qc_plots.tsv" % run)
    ret.append("analysis/report2/data_quality/coverage_statistics.tsv")
    return ret

rule report2_data_quality_order:
    """Dictate the ordering of these parts"""
    input:
        report2_data_quality_orderInputFn
    params:
        order= yaml_dump({'order': ['mapping_stats.tsv', '%s-qc_plots.tsv' % list(config['runs'].keys())[0], "coverage_statistics.tsv"]})
    output:
        "analysis/report2/data_quality/data_quality.yaml"
    group: "report2"
    shell:
        """echo '{params.order}' > {output}"""


###############################################################################

def report2_somatic_variants_summary_tblInputFn(wildcards):
    ls = []
    caller = config['somatic_caller']
    run = list(config['runs'].keys())[0]
    ls.append("analysis/somatic/somatic_mutation_summaries.%s.csv" % caller)
    ls.append("analysis/somatic/somatic_functional_annot_summaries.%s.csv" % caller)
    ls.append("analysis/somatic/%s/%s_%s_somatic_SNV_summaries.csv" % (run, run, caller))
    return ls

rule report2_somatic_variants_summary_tbls:
    input:
        report2_somatic_variants_summary_tblInputFn
    output:
        csv1 = "analysis/report2/somatic_variants/summary_table.csv",
        csv2 = "analysis/report2/somatic_variants/functional_annotation.csv",
        csv3 = "analysis/report2/somatic_variants/SNV_statistics.csv",
    shell:
        "cp {input[0]} {output.csv1} && cp {input[1]} {output.csv2} && cp {input[2]} {output.csv3}"

def report_legoPlotInputFn(wildcards):
    run = wildcards.run
    caller = config['somatic_caller']
    return "analysis/somatic/%s/%s_%s.filter.pdf" % (run, run, caller)

rule report2_somatic_variants_legoPlot:
    """Add the lego plot to the somatic_variants section of the report"""
    input:
        report_legoPlotInputFn
    output:
        png = "analysis/report2/somatic_variants/{run}_lego_plot.png",
        #sub_cap = "analysis/report2/somatic_variants/{run}_lego_plot_subcaption.txt",
        #cap = "analysis/report2/somatic_variants/{run}_lego_plot_caption.txt",
    params:
        page = 1,
        #sub_caption = "This is a test caption"
    message:
        "REPORT: generating somatic_variants lego plot"
    group: "report"
    shell:
        #"""echo "{params.sub_caption}" > {output.sub_cap} &&
        """Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output.png} {params.page}"""

def report2_somatic_variants_germlineCompareInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/germline/%s/%s_vcfcompare.txt" % (run, run)

rule report2_somatic_variants_germlineCompare:
    """report germline overlap"""
    input:
        report2_somatic_variants_germlineCompareInputFn
    params:
        run = list(config['runs'].keys())[0],
        sub = "NOTE: the % overlap was calculated using the Tumor TMB as the denominator",
    output:
        tsv = "analysis/report2/somatic_variants/tumor_mutational_burden.tsv",
        subcap = "analysis/report2/somatic_variants/tumor_mutational_burden_subcaption.txt",
    shell:
        """echo "{params.sub}" > {output.subcap} && cidc_wes/modules/scripts/report_somatic_tmb.py -f {input} -r {params.run} -o {output.tsv}"""

rule report2_somatic_variants_order:
    """Dictate the ordering of these parts"""
    input:
        "analysis/report2/somatic_variants/summary_table.csv",
        "analysis/report2/somatic_variants/functional_annotation.csv",
        "analysis/report2/somatic_variants/SNV_statistics.csv",
        #Lego plot
    params:
        order= yaml_dump({'order': ['summary_table.csv', 'functional_annotation.csv', 'SNV_statistics.csv', 'tumor_mutational_burden.tsv', "%s_lego_plot.png" % list(config['runs'].keys())[0]]})
    output:
        "analysis/report2/somatic_variants/somatic_variants.yaml"
    group: "report2"
    shell:
        """echo '{params.order}' > {output}"""
###############################################################################
def report2_data_quality_purityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run)

rule report2_data_quality_purity:
    """report tumor purity"""
    input:
        report2_data_quality_purityInputFn
    params:
        run = list(config['runs'].keys())[0]
    output:
        "analysis/report2/copy_number/tumor_purity.tsv"
    shell:
        """cidc_wes/modules/scripts/report_cnv_purity.py -f {input} -r {params.run} -o {output}"""

def report2_data_quality_clonalityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/clonality/%s/%s_table.tsv" % (run,run)

rule report2_data_quality_clonality:
    """report tumor clonality"""
    input:
        report2_data_quality_clonalityInputFn
    params:
        run = list(config['runs'].keys())[0]
    output:
        "analysis/report2/copy_number/tumor_clonality.tsv"
    shell:
        """cidc_wes/modules/scripts/report_cnv_clonality.py -f {input} -r {params.run} -o {output}"""

###############################################################################
def report2_neoantigens_HLAInputFn(wildcards):
    runName = list(config['runs'].keys())[0]
    run = config['runs'][runName]
    normal = run[0]
    tumor = run[1]
    ls = []
    if config.get('neoantigen_run_classII'):
        #optitype and xhla results
        ls = ["analysis/optitype/%s/%s_result.tsv" % (normal, normal),
               "analysis/xhla/%s/report-%s-hla.json" % (normal, normal),
              "analysis/optitype/%s/%s_result.tsv" % (tumor, tumor),
               "analysis/xhla/%s/report-%s-hla.json" % (tumor, tumor)]
    else:
        #optitype only
        ls = ["analysis/optitype/%s/%s_result.tsv" % (normal, normal),
              "analysis/optitype/%s/%s_result.tsv" % (tumor, tumor)]
    return ls

rule report2_neoantigens_HLA:
    """report HLA type"""
    input:
        report2_neoantigens_HLAInputFn
    params:
        normal = lambda wildcards, input: ",".join(input[:2]) if len(input) > 2 else input[0],
        tumor = lambda wildcards, input: ",".join(input[2:4]) if len(input) > 2 else input[1],
        names = ",".join(config['runs'][list(config['runs'].keys())[0]]),
    output:
        "analysis/report2/neoantigens/HLA_results.tsv"
    shell:
        """cidc_wes/modules/scripts/report_neoantigens_hla.py -n {params.normal} -t {params.tumor} -s {params.names} -o {output}"""

def report2_neoantigens_neoantigen_listInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/neoantigen/%s/%s_neoantigen_table.tsv" % (run, run)

rule report2_neoantigens_neoantigen_list:
    """report HLA type"""
    input:
        report2_neoantigens_neoantigen_listInputFn
    output:
        "analysis/report2/neoantigens/neoantigen_list.tsv"
    shell:
        "cut -f 1,3,4,5,6,7,8,8,10 {input} > {output}"

###############################################################################
rule report2_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report2_targets
    params:
        report_path = "analysis/report2",
        sections_list=",".join(['wes_meta','data_quality','somatic_variants','copy_number','neoantigens'])
    output:
        "analysis/report2/report.html"
    message:
        "REPORT: Generating WES report"
    group: "report2"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -o {output} && cp -r cidc_wes/report2/static {params.report_path}"""

