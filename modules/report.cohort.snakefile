# Len Taing 2020 (TGBTG)
#MODULE: wes cohort report module

import os
from yaml import dump as yaml_dump
from collections import OrderedDict

configfile: "config.cohort.yaml"

def cohort_report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #Meta information
    ls.append("analysis/cohort_report/runs_meta.json")
    ls.append("analysis/cohort_report/samples_meta.json")

    #Data Quality
    ls.append("analysis/cohort_report/data_quality/01_mapping_plots_bar.plotly")
    ls.append("analysis/cohort_report/data_quality/02_coverage_table.mqc")
    #ls.append("analysis/cohort_report/data_quality/02_mapping_stats.csv")
    ls.append("analysis/cohort_report/data_quality/03_gc_content_line.mqc")
    ls.append("analysis/cohort_report/data_quality/04_insert_size_line.mqc")

    #Copynumber
    ls.append("analysis/cohort_report/copy_number/01_copy_number_table.mqc")

    #Somatic
    ls.append("analysis/cohort_report/somatic/01_somatic_summary_table.mqc")
    ls.append("analysis/cohort_report/somatic/somatic_summary.json")
    ls.append("analysis/cohort_report/somatic/ti_tv.json")
    ls.append("analysis/cohort_report/somatic/tmb.json")
    ls.append("analysis/cohort_report/somatic/functional_summary.json")

    #Neoantigen
    ls.append("analysis/cohort_report/neoantigen/01_HLA_table.mqc")

    return ls

rule cohort_report_all:
    input:
        "analysis/cohort_report/report.html"

def cohort_report_inputFn(wildcards):
    """Returns a list of all of the json files"""
    ls = []
    for r in config['runs']:
        ls.append(config['runs'][r][0])
    return ls

###############################################################################
rule cohort_report_meta:
    """Generate the meta information for the report"""
    input:
        cohort_report_inputFn
    output:
        runs="analysis/cohort_report/runs_meta.json",
        samples="analysis/cohort_report/samples_meta.json",
    params:
        files = lambda wildcards,input: " -f ".join(input),
    message:
        "REPORT: creating meta json files"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/cohort_report/cr_meta.py -f {params.files} -r {output.runs} -s {output.samples}"""

###############################################################################
rule cohort_report_data_quality_plots:
    """Generate the mapping stats plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/01_mapping_plots_bar.plotly",
        details="analysis/cohort_report/data_quality/01_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This table shows the total number reads in each sample, how many of those reads were mapped, and how many are de-duplicated reads.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0}}),
    message:
        "REPORT: creating mapping plots for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_mappingPlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_coverage_table:
    """Generate the coverage table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/02_coverage_table.mqc",
        details="analysis/cohort_report/data_quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This table shows read depth coverage of each sample.'""",
        plot_options = "cpswitch: False",
    message:
        "REPORT: creating coverage table for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_coverageTable.py -f {params.files} -o {output.csv}"""

# rule cohort_report_data_quality_table:
#     """Generate the mapping stats table for the report"""
#     input:
#         cohort_report_inputFn
#     output:
#         csv="analysis/cohort_report/data_quality/02_mapping_stats.csv",
#         details="analysis/cohort_report/data_quality/02_details.yaml",
#     params:
#         files = lambda wildcards,input: " -f ".join(input),
#         caption="""caption: 'This table shows the total number reads in each sample, how many of those reads were mapped, and how many are de-duplicated reads.'"""
#     message:
#         "REPORT: creating mapping stats for data_quality section"
#     group: "cohort_report"
#     shell:
#         """echo "{params.caption}" >> {output.details} && cidc_wes/modules/scripts/cohort_report/cr_dataQual_mappingStats.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_gc_plots:
    """Generate the gc content plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/03_gc_content_line.mqc",
        details="analysis/cohort_report/data_quality/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        subcaption="""subcaption: 'GC Plot shows the distribution of %GC bases within a 100bp window. In human, the mean GC content is approx. 40%.'""",
    message:
        "REPORT: creating gc content plots for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.subcaption}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_gcPlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_insertSize_plots:
    """Generate the insert size plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/04_insert_size_line.mqc",
        #details="analysis/cohort_report/data_quality/04_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        #subcaption="""subcaption: 'GC Plot shows the distribution of %GC bases within a 100bp window. In human, the mean GC content is approx. 40%.'""",
    message:
        "REPORT: creating insert size plots for data_quality section"
    group: "cohort_report"
    shell:
        #"""echo "{params.subcaption}" >> {output.details} && 
        """cidc_wes/modules/scripts/cohort_report/cr_dataQual_insertSizePlots.py -f {params.files} -o {output.csv}"""
###############################################################################
rule cohort_report_copynumber_table:
    """Generate the copynumber table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/copy_number/01_copy_number_table.mqc",
        #details="analysis/cohort_report/data_quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        #caption="""caption: 'This table shows read depth coverage of each sample.'""",
        #plot_options = "cpswitch: False",
    message:
        "REPORT: creating copynumber table for copy_number section"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -o {output.csv}"""
###############################################################################
rule cohort_report_somatic_summary_table:
    """Generate the somatic summary table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/somatic/01_somatic_summary_table.mqc",
        ss="analysis/cohort_report/somatic/somatic_summary.json",
        titv="analysis/cohort_report/somatic/ti_tv.json",
        tmb="analysis/cohort_report/somatic/tmb.json",
        func="analysis/cohort_report/somatic/functional_summary.json",
        #...other json here
        #details="analysis/cohort_report/data_quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        #caption="""caption: 'This table shows read depth coverage of each sample.'""",
        #plot_options = "cpswitch: False",
    message:
        "REPORT: creating somatic summary table for somatic section"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/cohort_report/cr_somatic_somaticSummaryTable.py -f {params.files} -o {output.csv} -j {output.ss} -k {output.titv} -l {output.tmb} -m {output.func}"""

###############################################################################
def getHLATable_categories():
    #Actual
    #classII = ['DPB1-1', 'DPB1-2', 'DQB1-1', 'DQB1-2', 'DRB1-1', 'DRB1-2']
    #Prettyprint version
    classI = list(map(lambda x: x.title(), ['A-1', 'A-2', 'B-1', 'B-2', 'C-1', 'C-2']))
    classII = list(map(lambda x: x.title(), ['DPB1-1', 'DPB1-2', 'DQB1-1', 'DQB1-2', 'DRB1-1', 'DRB1-2']))
    both = classI + classII
    tmp = {}
    for c in both:
        if c in classI:
            tmp[c] = {'hidden': False}
        else:
            tmp[c] = {'hidden': True}
    ret = yaml_dump({'cats': tmp })
    #print(ret)
    return ret

rule cohort_report_HLA_table:
    """Generate the HLA table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/neoantigen/01_HLA_table.mqc",
        details="analysis/cohort_report/neoantigen/01_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        #caption="""caption: 'Table of HLA Alleles.\n**NOTE: Click on Configure Table to show Class II alleles**'""",
        caption="""caption: '**NOTE: Click on Configure Table to show Class II alleles**'""",
        table_options = "table_title: 'HLA Alleles Table'",
        #Turn off display of Class I Alleles
        cats = getHLATable_categories(),
    message:
        "REPORT: creating HLA table for neoantigen section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        echo "{params.table_options}" >> {output.details} &&
        echo "{params.cats}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_neoantigen_hlaTable.py -f {params.files} -o {output.csv}"""

###############################################################################
rule cohort_report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        cohort_report_targets
    params:
        jinja2_template="cidc_wes/report2/index.cohort.html",
        report_path = "analysis/cohort_report",
        sections_list=",".join(['data_quality','somatic', 'copy_number', 'neoantigen'])
        #sections_list=",".join(['somatic'])
    output:
        "analysis/cohort_report/report.html"
    message:
        "REPORT: Generating WES report"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -t {params.jinja2_template} -o {output} && cp -r cidc_wes/report2/static {params.report_path}"""

