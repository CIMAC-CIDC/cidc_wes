# Len Taing 2020 (TGBTG)
#MODULE: wes cohort report module

import os
from yaml import dump as yaml_dump

configfile: "config.cohort.yaml"

def cohort_report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #Data Quality
    ls.append("analysis/cohort_report/data_quality/01_mapping_plots_table.plot")
    ls.append("analysis/cohort_report/data_quality/02_coverage_table.plot")
    #ls.append("analysis/cohort_report/data_quality/02_mapping_stats.csv")
    ls.append("analysis/cohort_report/data_quality/03_gc_content_line.plot")
    ls.append("analysis/cohort_report/data_quality/04_insert_size_line.plot")

    #Copynumber
    ls.append("analysis/cohort_report/copy_number/01_copy_number_table.plot")
    
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
rule cohort_report_data_quality_plots:
    """Generate the mapping stats plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/01_mapping_plots_table.plot",
        details="analysis/cohort_report/data_quality/01_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This table shows the total number reads in each sample, how many of those reads were mapped, and how many are de-duplicated reads.'""",
        plot_options = "cpswitch: False",
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
        csv="analysis/cohort_report/data_quality/02_coverage_table.plot",
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
        csv="analysis/cohort_report/data_quality/03_gc_content_line.plot",
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
        csv="analysis/cohort_report/data_quality/04_insert_size_line.plot",
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
        csv="analysis/cohort_report/copy_number/01_copy_number_table.plot",
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
rule cohort_report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        cohort_report_targets
    params:
        report_path = "analysis/cohort_report",
        sections_list=",".join(['data_quality','copy_number'])
    output:
        "analysis/cohort_report/report.html"
    message:
        "REPORT: Generating WES report"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -o {output} && cp -r cidc_wes/report2/static {params.report_path}"""

