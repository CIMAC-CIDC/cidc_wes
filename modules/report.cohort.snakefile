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
    ls.append("analysis/cohort_report/data_quality/01_mapping_bar.plotly")
    ls.append("analysis/cohort_report/data_quality/02_coverage_bar.plotly")
    ls.append("analysis/cohort_report/data_quality/03_gc_content_line.plotly")
    ls.append("analysis/cohort_report/data_quality/04_insert_size_line.plotly")
    ls.append("analysis/cohort_report/data_quality/05_mean_quality_bar.plotly")

    #Copynumber
    ls.append("analysis/cohort_report/copy_number/01_clonality_bar.plotly")
    ls.append("analysis/cohort_report/copy_number/02_purity_bar.plotly")
    ls.append("analysis/cohort_report/copy_number/03_ploidy_line.plotly")

    #Somatic
    ls.append("analysis/cohort_report/somatic/cohort.maf.gz")
    ls.append("analysis/cohort_report/somatic/01_somatic_variants_summary.png")
    ls.append("analysis/cohort_report/somatic/02_oncoplot.png")
    ls.append("analysis/cohort_report/somatic/03_Ti-Tv.png")
    ls.append("analysis/cohort_report/somatic/04_vaf.png")
    ls.append("analysis/cohort_report/somatic/05_tcga_comparison.png")
    ls.append("analysis/cohort_report/somatic/06_somatic_interactions.png")
    ls.append("analysis/cohort_report/somatic/07_variants_oncoplot_oncoplot.plotly")
    ls.append("analysis/cohort_report/somatic/10_lollipop_plots.csv")
    
    ls.append("analysis/cohort_report/somatic/plots/20_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic/plots/21_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic/plots/22_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic/plots/23_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic/plots/24_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic/plots/25_lollipop_plot.png")

    # ls.append("analysis/cohort_report/somatic/01_mutation_summary_bar.plotly")
    # ls.append("analysis/cohort_report/somatic/02_somatic_summary_table.mqc")
    # ls.append("analysis/cohort_report/somatic/somatic_summary.json")
    # ls.append("analysis/cohort_report/somatic/ti_tv.json")
    # ls.append("analysis/cohort_report/somatic/tmb.json")
    # ls.append("analysis/cohort_report/somatic/functional_summary.json")

    #HLA
    ls.append("analysis/cohort_report/HLA/01_HLA_table.dt")
    ls.append("analysis/cohort_report/HLA/02_HLA_histogram_histogram.plotly")
    ls.append("analysis/cohort_report/HLA/03_HLA_oncoplot_oncoplot.plotly")
    #ls.append("analysis/cohort_report/neoantigen/03_neoantigen_table.csv")
    #ls.append("analysis/cohort_report/neoantigen/neoantigen_table.json")
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
        csv="analysis/cohort_report/data_quality/01_mapping_bar.plotly",
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
        csv="analysis/cohort_report/data_quality/02_coverage_bar.plotly",
        details="analysis/cohort_report/data_quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This table shows read depth coverage of each sample.'""",
        plot_options = yaml_dump({'plotly': {'x':'Percent Bases >50', 'hover_data':['Q1 Depth','Mean Depth','Median Depth','Q3 Depth']}}),
    message:
        "REPORT: creating coverage table for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_coverageTable.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_gc_plots:
    """Generate the gc content plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/03_gc_content_line.plotly",
        details="analysis/cohort_report/data_quality/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        subcaption="""subcaption: 'GC Plot shows the distribution of %GC bases within a 100bp window. In human, the mean GC content is approx. 40%.'""",
        plot_options = yaml_dump({'plotly': {'labels':{'X':'% GC bases','value':''}}}),
    message:
        "REPORT: creating gc content plots for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.subcaption}" >> {output.details} && 
        echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_gcPlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_insertSize_plots:
    """Generate the insert size plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/04_insert_size_line.plotly",
        details="analysis/cohort_report/data_quality/04_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        plot_options = yaml_dump({'plotly': {'labels':{'X':'Insert Size','value':'Counts'}}}),
    message:
        "REPORT: creating insert size plots for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.plot_options}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_insertSizePlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_meanQual:
    """Generate the mean quality plots for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/data_quality/05_mean_quality_bar.plotly",
        details="analysis/cohort_report/data_quality/05_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        plot_options = yaml_dump({'plotly': {'color_discrete_sequence':["#44bd32"]}}), #make the bars green
    message:
        "REPORT: creating mean quality for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_dataQual_meanQual.py -f {params.files} -o {output.csv}"""

###############################################################################
rule cohort_report_copynumber_clonality:
    """Generate the clonality table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/copy_number/01_clonality_bar.plotly",
        #details="analysis/cohort_report/data_quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        attr="clonality",
        #caption="""caption: 'This table shows read depth coverage of each sample.'""",
        #plot_options = "cpswitch: False",
    message:
        "REPORT: creating clonality plot for copy_number section"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -a {params.attr} -o {output.csv}"""

rule cohort_report_copynumber_purity:
    """Generate the purity plot for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/copy_number/02_purity_bar.plotly",
        details="analysis/cohort_report/copy_number/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        attr="purity",
        plot_options = yaml_dump({'plotly': {'color_discrete_sequence':["#0097e6"]}}), #make the bars blue
    message:
        "REPORT: creating copynumber table for copy_number section"
    group: "cohort_report"
    shell:
        """echo "{params.plot_options}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -a {params.attr} -o {output.csv}"""

rule cohort_report_copynumber_ploidy:
    """Generate the ploidy plot for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/copy_number/03_ploidy_line.plotly",
        details="analysis/cohort_report/copy_number/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        attr=" -a ".join(["ploidy","dipLogR"]),
        plot_options = yaml_dump({'plotly': {'color_discrete_sequence':["#e1b12c",'#7f8fa6']}}), #make the lines yellow and gray
    message:
        "REPORT: creating copynumber table for copy_number section"
    group: "cohort_report"
    shell:
        """echo "{params.plot_options}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -a {params.attr} -o {output.csv}"""

###############################################################################
rule cohort_report_somatic_makeCohort:
    """Collect the maf files and gzip them for mafTools input"""
    input:
        cohort_report_inputFn
    output:
        #LATER: make this a temp file
        gz=temp("analysis/cohort_report/somatic/cohort.maf.gz"),
    params:
        files = lambda wildcards,input: " -f ".join(input),
    message:
        "REPORT: creating cohort maf bundle"
    group: "cohort_report"
    shell:
        #NOTE: this code needs speedup!
        """cidc_wes/modules/scripts/cohort_report/cr_somatic_bundleCohortMafs.py -f {params.files} -o {output.gz}"""

rule cohort_report_somatic_mafTools:
    """Generate the mafTools plots for the report"""
    input:
        mafs="analysis/cohort_report/somatic/cohort.maf.gz",
        cancerGeneList = "cidc_wes/static/oncoKB/cancerGeneList.tsv",
    output:
        summary="analysis/cohort_report/somatic/01_somatic_variants_summary.png",
        onco="analysis/cohort_report/somatic/02_oncoplot.png",
        titv="analysis/cohort_report/somatic/03_Ti-Tv.png",
        vaf="analysis/cohort_report/somatic/04_vaf.png",
        tcga="analysis/cohort_report/somatic/05_tcga_comparison.png",
        interact="analysis/cohort_report/somatic/06_somatic_interactions.png",

        lolli01= "analysis/cohort_report/somatic/plots/20_lollipop_plot.png",
        lolli02= "analysis/cohort_report/somatic/plots/21_lollipop_plot.png",
        lolli03= "analysis/cohort_report/somatic/plots/22_lollipop_plot.png",
        lolli04= "analysis/cohort_report/somatic/plots/23_lollipop_plot.png",
        lolli05= "analysis/cohort_report/somatic/plots/24_lollipop_plot.png",
        lolli06= "analysis/cohort_report/somatic/plots/25_lollipop_plot.png",

        #details="analysis/cohort_report/somatic/01_details.yaml",
    params:
        #caption="""caption: 'This table shows read depth coverage of each sample.'""",
    message:
        "REPORT: creating mafTools plot for somatic section"
    group: "cohort_report"
    shell:
        """Rscript cidc_wes/modules/scripts/cohort_report/cr_somatic_mafPlots.R {input.mafs} {input.cancerGeneList} {output.summary} {output.onco} {output.titv} {output.vaf} {output.tcga} {output.interact} {output.lolli01} {output.lolli02} {output.lolli03} {output.lolli04} {output.lolli05} {output.lolli06}"""

rule cohort_report_somatic_lollipop_table:
    """Generate the table of lollipop plots"""
    input:
        lolli01= "analysis/cohort_report/somatic/plots/20_lollipop_plot.png",
        lolli02= "analysis/cohort_report/somatic/plots/21_lollipop_plot.png",
        lolli03= "analysis/cohort_report/somatic/plots/22_lollipop_plot.png",
        lolli04= "analysis/cohort_report/somatic/plots/23_lollipop_plot.png",
        lolli05= "analysis/cohort_report/somatic/plots/24_lollipop_plot.png",
        lolli06= "analysis/cohort_report/somatic/plots/25_lollipop_plot.png",
    output:
        csv = "analysis/cohort_report/somatic/10_lollipop_plots.csv",
        details="analysis/cohort_report/somatic/10_details.yaml",
    params:
        caption="""caption: 'This table shows lollipop plots for the top 6 cancer driving genes.  The first row shows the top 3 genes while the second row shows genes 4 through 6.'""",
        report_path = "analysis/cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&      
        cidc_wes/modules/scripts/cohort_report/cr_somatic_lollipop_table.py -f {input.lolli01} -f {input.lolli02} -f {input.lolli03} -f {input.lolli04} -f {input.lolli05} -f {input.lolli06} -r {params.report_path} -o {output}"""
    
rule cohort_report_somatic_variants_oncoplot:
    """Generate the variants based oncoplot for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/somatic/07_variants_oncoplot_oncoplot.plotly",
        details="analysis/cohort_report/somatic/07_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'Oncoplot of shared somatic variants'""",
        plot_options = yaml_dump({'plotly': {'top_ngenes':25, 'colors': ['#b5b5b5', '#8c7ae6']}}),
    message:
        "REPORT: creating variants based oncoplot for somatic section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
          echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_somatic_variantOncoplot.py -f {params.files} -o {output.csv}"""

###############################################################################
rule cohort_report_HLA_table:
    """Generate the HLA table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/HLA/01_HLA_table.dt",
        details="analysis/cohort_report/HLA/01_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        table_options = "table_title: 'HLA Alleles Table'",
    message:
        "REPORT: creating HLA table for HLA section"
    group: "cohort_report"
    shell:
        """echo "{params.table_options}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_hla_hlaTable.py -f {params.files} -o {output.csv}"""

rule cohort_report_HLA_histogram:
    """Generate the HLA histogram plot for the report"""
    input:
        csv="analysis/cohort_report/HLA/01_HLA_table.dt",
    output:
        csv="analysis/cohort_report/HLA/02_HLA_histogram_histogram.plotly",
        details="analysis/cohort_report/HLA/02_details.yaml",
    params:
        caption="""caption: 'Histogram of shared HLA-alleles'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0}}),
    message:
        "REPORT: creating HLA table for HLA section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        echo "{params.plot_options}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_hla_hlaHistogram.py -f {input} -o {output.csv}"""

rule cohort_report_HLA_oncoplot:
    """Generate the HLA histogram plot for the report"""
    input:
        csv="analysis/cohort_report/HLA/01_HLA_table.dt",
    output:
        csv="analysis/cohort_report/HLA/03_HLA_oncoplot_oncoplot.plotly",
        details="analysis/cohort_report/HLA/03_details.yaml",
    params:
        caption="""caption: 'Oncoplot of shared HLA-alleles'""",
        plot_options = yaml_dump({'plotly': {'top_ngenes':15, 'colors': ['#b5b5b5', '#e84118']}}),
    message:
        "REPORT: creating HLA oncoplot for HLA section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
          echo "{params.plot_options}" >> {output.details} && 
        cidc_wes/modules/scripts/cohort_report/cr_hla_hlaOncoplot.py -f {input} -o {output.csv}"""

#DEPRECATED
rule cohort_report_neoantigen_table:
    """Generate the neoantigen table for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/neoantigen/03_neoantigen_table.csv",
        json="analysis/cohort_report/neoantigen/neoantigen_table.json",
        details="analysis/cohort_report/neoantigen/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'Table of neoantigens for each sample.  Only the top hit is shown.  Select a row and click the download button to save the full list to disk.'""",
    message:
        "REPORT: creating neonatigen table for neoantigen section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_neoantigen_neoantigenTable.py -f {params.files} -o {output.csv} -j {output.json}"""

###############################################################################
rule cohort_report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        cohort_report_targets
    params:
        jinja2_template="cidc_wes/report/index.cohort.html",
        report_path = "analysis/cohort_report",
        sections_list=",".join(['data_quality','copy_number', 'somatic', 'HLA']),
        title="WES Cohort Report",
        #sections_list=",".join(['somatic'])
    output:
        "analysis/cohort_report/report.html"
    message:
        "REPORT: Generating WES report"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -j {params.jinja2_template} -t "{params.title}" -o {output} && cp -r cidc_wes/report/static {params.report_path}"""

