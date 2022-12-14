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
    #ls.append("analysis/cohort_report/samples_meta.json")
    ls.append("analysis/cohort_report/wes_data.json")

    #Data Quality
    ls.append("analysis/cohort_report/data_quality/01_mapping.stub")
    ls.append("analysis/cohort_report/data_quality/02_coverage.stub")
    ls.append("analysis/cohort_report/data_quality/03_GC_Content.stub")
    ls.append("analysis/cohort_report/data_quality/04_insert_size.stub")
    ls.append("analysis/cohort_report/data_quality/05_mean_quality.stub")

    #Copynumber
    ls.append("analysis/cohort_report/copy_number_variation/01_clonality.stub")
    ls.append("analysis/cohort_report/copy_number_variation/02_purity.stub")
    ls.append("analysis/cohort_report/copy_number_variation/03_ploidy.stub")

    #Somatic
    ls.append("analysis/cohort_report/somatic_variants/cohort.maf.gz")
    #summary plots
    ls.append("analysis/cohort_report/somatic_variants/01_variants_per_sample.stub")
    ls.append("analysis/cohort_report/somatic_variants/02_tumor_mutational_burden.stub")
    ls.append("analysis/cohort_report/somatic_variants/03_top_10_genes.stub")
    ls.append("analysis/cohort_report/somatic_variants/04_variant_class.stub")
    ls.append("analysis/cohort_report/somatic_variants/05_variant_class_summary.stub")
    ls.append("analysis/cohort_report/somatic_variants/06_variant_type.stub")
    ls.append("analysis/cohort_report/somatic_variants/07_SNV_Class.stub")
    
    ls.append("analysis/cohort_report/somatic_variants/08_Ti-Tv.stub")
    ls.append("analysis/cohort_report/somatic_variants/09_lego_plot.stub")
    ls.append("analysis/cohort_report/somatic_variants/10_oncoplot.stub")
    ls.append("analysis/cohort_report/somatic_variants/11_VAF.png")
    ls.append("analysis/cohort_report/somatic_variants/12_TCGA_Comparison.png")
    ls.append("analysis/cohort_report/somatic_variants/13_somatic_interactions.png")
    #Removing this plot b/c it's looking strange
    #ls.append("analysis/cohort_report/somatic_variants/07_variants_oncoplot_oncoplot.plotly")
    ls.append("analysis/cohort_report/somatic_variants/14_lollipop_plots.csv")
    
    ls.append("analysis/cohort_report/somatic_variants/plots/20_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic_variants/plots/21_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic_variants/plots/22_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic_variants/plots/23_lollipop_plot.png")
    ls.append("analysis/cohort_report/somatic_variants/plots/24_lollipop_plot.png")

    # ls.append("analysis/cohort_report/somatic/01_mutation_summary_bar.plotly")
    # ls.append("analysis/cohort_report/somatic/02_somatic_summary_table.mqc")
    # ls.append("analysis/cohort_report/somatic/somatic_summary.json")
    # ls.append("analysis/cohort_report/somatic/ti_tv.json")
    # ls.append("analysis/cohort_report/somatic/tmb.json")
    # ls.append("analysis/cohort_report/somatic/functional_summary.json")

    #HLA
    ls.append("analysis/cohort_report/HLA/01_HLA_Oncoplot.stub")
    ls.append("analysis/cohort_report/HLA/02_HLA_Table.dt")
    #ls.append("analysis/cohort_report/HLA/02_HLA_Histogram_histogram.plotly")

    #ls.append("analysis/cohort_report/neoantigen/03_neoantigen_table.csv")
    #ls.append("analysis/cohort_report/neoantigen/neoantigen_table.json")

    #MISC
    ls.append("analysis/cohort_report/MISC/01_tcellextrect.stub")
    ls.append("analysis/cohort_report/MISC/02_msisensor2.stub")    
    return ls

rule cohort_report_all:
    input:
        "analysis/cohort_report/report.html"
    benchmark: "benchmarks/cohort_report/cohort_report_all.txt"

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
        jsons=cohort_report_inputFn,
        metasheet="metasheet.cohort.csv",
    output:
        runs="analysis/cohort_report/runs_meta.json",
        #samples="analysis/cohort_report/samples_meta.json",
        wes_data="analysis/cohort_report/wes_data.json",
    params:
        files = lambda wildcards,input: " -f ".join(input.jsons),
    message:
        "REPORT: creating json files"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/cohort_report/cr_meta.py -f {params.files} -m {input.metasheet} -r {output.runs} -d {output.wes_data}"""

###############################################################################
rule cohort_report_data_quality_plots:
    """Generate the mapping stats plots for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/data_quality/01_mapping_bar.plotly",
        csv="analysis/cohort_report/data_quality/01_mapping.stub",
        details="analysis/cohort_report/data_quality/01_details.yaml",
    params:
        #files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This plot shows the total number reads in each sample, how many of those reads were mapped, and how many are de-duplicated reads.'""",
        #plot_options = yaml_dump({'render': False, 'plotly': {'barmode':"overlay",'opacity':1.0}}),
    message:
        "REPORT: creating mapping plots for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        touch {output.csv}
        """
        #cidc_wes/modules/scripts/cohort_report/cr_dataQual_mappingPlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_coverage_table:
    """Generate the coverage table for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/data_quality/02_coverage_bar.plotly",
        csv="analysis/cohort_report/data_quality/02_coverage.stub",
        details="analysis/cohort_report/data_quality/02_details.yaml",
    params:
        #files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This plot shows read depth coverage of each sample. Hover over each bar to see more coverage information for the associated sample.'""",
        #plot_options = yaml_dump({'render': False, 'plotly': {'x':'Mean Depth', 'hover_data':['Q1 Depth','Median Depth','Q3 Depth', 'Percent Bases >50']}}),
    message:
        "REPORT: creating coverage table for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        touch {output.csv}
        """
        #cidc_wes/modules/scripts/cohort_report/cr_dataQual_coverageTable.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_gc_plots:
    """Generate the gc content plots for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/data_quality/03_GC_Content_line.plotly",
        csv="analysis/cohort_report/data_quality/03_GC_Content.stub",
        details="analysis/cohort_report/data_quality/03_details.yaml",
    params:
        #files = lambda wildcards,input: " -f ".join(input),
        subcaption="""subcaption: 'GC Plot shows the distribution of %GC bases within a 100bp window. In human, the mean GC content is approx. 40%.'""",
        #plot_options = yaml_dump({'render': False, 'plotly': {'labels':{'X':'% GC bases','value':''}}}),
    message:
        "REPORT: creating gc content plots for data_quality section"
    group: "cohort_report"
    shell:
        """echo "{params.subcaption}" >> {output.details} && 
        touch {output.csv}"""
        #echo "{params.plot_options}" >> {output.details} && 
        #cidc_wes/modules/scripts/cohort_report/cr_dataQual_gcPlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_insertSize_plots:
    """Generate the insert size plots for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/data_quality/04_insert_size_line.plotly",
        csv="analysis/cohort_report/data_quality/04_insert_size.stub",
        #details="analysis/cohort_report/data_quality/04_details.yaml",
    #params:
        #files = lambda wildcards,input: " -f ".join(input),
        #plot_options = yaml_dump({'render': False, 'plotly': {'labels':{'X':'Insert Size','value':'Counts'}}}),
    message:
        "REPORT: creating insert size plots for data_quality section"
    group: "cohort_report"
    shell:
        """ touch {output.csv}
        """
        #echo "{params.plot_options}" >> {output.details} &&
        #cidc_wes/modules/scripts/cohort_report/cr_dataQual_insertSizePlots.py -f {params.files} -o {output.csv}"""

rule cohort_report_data_quality_meanQual:
    """Generate the mean quality plots for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/data_quality/05_mean_quality_bar.plotly",
        csv="analysis/cohort_report/data_quality/05_mean_quality.stub",
        #details="analysis/cohort_report/data_quality/05_details.yaml",
    params:
        #files = lambda wildcards,input: " -f ".join(input),
        #plot_options = yaml_dump({'render': False, 'plotly': {'color_discrete_sequence':["#44bd32"]}}), #make the bars green
    message:
        "REPORT: creating mean quality for data_quality section"
    group: "cohort_report"
    shell:
        """touch {output.csv}"""
        #echo "{params.plot_options}" >> {output.details} && 
        #cidc_wes/modules/scripts/cohort_report/cr_dataQual_meanQual.py -f {params.files} -o {output.csv}"""

###############################################################################
rule cohort_report_copynumber_clonality:
    """Generate the clonality table for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        csv="analysis/cohort_report/copy_number_variation/01_clonality.stub",
        details="analysis/cohort_report/copy_number_variation/01_details.yaml",
    params:
        #files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'This plot shows the cancer cell fractions (CCF) of each cluster for each sample.'""",
        subcaption="""subcaption: '<b>NOTE: These cluster estimates were derived from single-sample runs, therefore clusters do not corresponde across samples, i.e. the cluster0 in sampleX does not correspond to cluster0 in sampleY.  Also, the CCFs do not add up to 1.0.</b>'""",
        #plot_options = yaml_dump({'render': False}),
    message:
        "REPORT: creating clonality plot for copy_number_variation section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && echo "{params.subcaption}" >> {output.details} && touch {output.csv}"""
        #cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -a {params.attr} -o {output.csv}"""

rule cohort_report_copynumber_purity:
    """Generate the purity plot for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/copy_number_variation/02_purity_bar.plotly",
        csv="analysis/cohort_report/copy_number_variation/02_purity.stub",
        #details="analysis/cohort_report/copy_number_variation/02_details.yaml",
    #params:
        #files = lambda wildcards,input: " -f ".join(input),
        #attr="purity",
        #plot_options = yaml_dump({'render': False, 'plotly': {'color_discrete_sequence':["#0097e6"]}}), #make the bars blue
    message:
        "REPORT: creating copynumber table for copy_number_variation section"
    group: "cohort_report"
    shell:
        """touch {output.csv}"""
        #echo "{params.plot_options}" >> {output.details} &&
        #cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -a {params.attr} -o {output.csv}"""

rule cohort_report_copynumber_ploidy:
    """Generate the ploidy plot for the report"""
    #input:
    #    cohort_report_inputFn
    output:
        #csv="analysis/cohort_report/copy_number_variation/03_ploidy_line.plotly",
        csv="analysis/cohort_report/copy_number_variation/03_ploidy.stub",
        #details="analysis/cohort_report/copy_number_variation/03_details.yaml",
    #params:
        #files = lambda wildcards,input: " -f ".join(input),
        #attr="ploidy", #-a ".join(["ploidy"]),
        #plot_options = yaml_dump({'render': False, 'plotly': {'color_discrete_sequence':["#e1b12c",'#7f8fa6']}}), #make the lines yellow and gray
    message:
        "REPORT: creating copynumber table for copy_number_variation section"
    group: "cohort_report"
    shell:
        """touch {output.csv}"""
        #echo "{params.plot_options}" >> {output.details} &&
        #cidc_wes/modules/scripts/cohort_report/cr_copynumber_cnvTable.py -f {params.files} -a {params.attr} -o {output.csv}"""

###############################################################################
rule cohort_report_somatic_makeCohort:
    """Collect the maf files and gzip them for mafTools input"""
    input:
        cohort_report_inputFn
    output:
        #LATER: make this a temp file
        gz=temp("analysis/cohort_report/somatic_variants/cohort.maf.gz"),
    params:
        files = lambda wildcards,input: " -f ".join(input),
    message:
        "REPORT: creating cohort maf bundle"
    group: "cohort_report"
    shell:
        #NOTE: this code needs speedup!
        """cidc_wes/modules/scripts/cohort_report/cr_somatic_bundleCohortMafs.py -f {params.files} -o {output.gz}"""

rule cohort_report_somatic_variant_stubs:
    input:
        mafs="analysis/cohort_report/somatic_variants/cohort.maf.gz",
        cancerGeneList = "cidc_wes/static/oncoKB/cancerGeneList.tsv",
    output: 
        varPerSample="analysis/cohort_report/somatic_variants/01_variants_per_sample.stub",
        tmb="analysis/cohort_report/somatic_variants/02_tumor_mutational_burden.stub",
        top10="analysis/cohort_report/somatic_variants/03_top_10_genes.stub",
        varClass="analysis/cohort_report/somatic_variants/04_variant_class.stub",
        varClassSum="analysis/cohort_report/somatic_variants/05_variant_class_summary.stub",
        varType="analysis/cohort_report/somatic_variants/06_variant_type.stub",
        snvClass="analysis/cohort_report/somatic_variants/07_SNV_Class.stub",
        titv="analysis/cohort_report/somatic_variants/08_Ti-Tv.stub",
        lego="analysis/cohort_report/somatic_variants/09_lego_plot.stub",
        onco="analysis/cohort_report/somatic_variants/10_oncoplot.stub",
    message:
        "REPORT: creating mafTools plot for somatic section"
    group: "cohort_report"
    shell:
        """touch {output.varPerSample} && touch {output.tmb} &&
        touch {output.top10} &&
        touch {output.varClass} && touch {output.varClassSum} &&
        touch {output.varType} && touch {output.snvClass} && 
        touch {output.titv} && touch {output.lego} && touch {output.onco}
        """
    
rule cohort_report_somatic_mafTools:
    """Generate the mafTools plots for the report"""
    input:
        mafs="analysis/cohort_report/somatic_variants/cohort.maf.gz",
        cancerGeneList = "cidc_wes/static/oncoKB/cancerGeneList.tsv",
    output:
        vaf="analysis/cohort_report/somatic_variants/11_VAF.png",
        tcga="analysis/cohort_report/somatic_variants/12_TCGA_Comparison.png",
        interact="analysis/cohort_report/somatic_variants/13_somatic_interactions.png",

        lolli01= "analysis/cohort_report/somatic_variants/plots/20_lollipop_plot.png",
        lolli02= "analysis/cohort_report/somatic_variants/plots/21_lollipop_plot.png",
        lolli03= "analysis/cohort_report/somatic_variants/plots/22_lollipop_plot.png",
        lolli04= "analysis/cohort_report/somatic_variants/plots/23_lollipop_plot.png",
        lolli05= "analysis/cohort_report/somatic_variants/plots/24_lollipop_plot.png",

        #details="analysis/cohort_report/somatic_variants/01_details.yaml",
    params:
        #caption="""caption: 'This table shows read depth coverage of each sample.'""",
    message:
        "REPORT: creating mafTools plot for somatic section"
    group: "cohort_report"
    shell:
        """Rscript cidc_wes/modules/scripts/cohort_report/cr_somatic_mafPlots.R {input.mafs} {input.cancerGeneList} {output.vaf} {output.tcga} {output.interact} {output.lolli01} {output.lolli02} {output.lolli03} {output.lolli04} {output.lolli05}"""

rule cohort_report_somatic_lollipop_table:
    """Generate the table of lollipop plots"""
    input:
        lolli01= "analysis/cohort_report/somatic_variants/plots/20_lollipop_plot.png",
        lolli02= "analysis/cohort_report/somatic_variants/plots/21_lollipop_plot.png",
        lolli03= "analysis/cohort_report/somatic_variants/plots/22_lollipop_plot.png",
        lolli04= "analysis/cohort_report/somatic_variants/plots/23_lollipop_plot.png",
        lolli05= "analysis/cohort_report/somatic_variants/plots/24_lollipop_plot.png",
    output:
        csv = "analysis/cohort_report/somatic_variants/14_lollipop_plots.csv",
        details="analysis/cohort_report/somatic_variants/14_details.yaml",
    params:
        caption="""caption: 'This table shows lollipop plots for the top 5 cancer driving genes.'""",
        report_path = "analysis/cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&      
        cidc_wes/modules/scripts/cohort_report/cr_somatic_lollipop_table.py -f {input.lolli01} -f {input.lolli02} -f {input.lolli03} -f {input.lolli04} -f {input.lolli05} -r {params.report_path} -o {output}"""

#THIS RULE is DEPRECATED
rule cohort_report_somatic_variants_oncoplot:
    """Generate the variants based oncoplot for the report"""
    input:
        cohort_report_inputFn
    output:
        csv="analysis/cohort_report/somatic_variants/07_variants_oncoplot_oncoplot.plotly",
        details="analysis/cohort_report/somatic_variants/07_details.yaml",
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
        csv="analysis/cohort_report/HLA/02_HLA_Table.dt",
        details="analysis/cohort_report/HLA/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        table_options = "table_title: 'HLA Alleles Table'",
    message:
        "REPORT: creating HLA table for HLA section"
    group: "cohort_report"
    shell:
        """echo "{params.table_options}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_hla_hlaTable.py -f {params.files} -o {output.csv}"""

#REMOVED!
# rule cohort_report_HLA_histogram:
#     """Generate the HLA histogram plot for the report"""
#     input:
#         csv="analysis/cohort_report/HLA/01_HLA_Table.dt",
#     output:
#         csv="analysis/cohort_report/HLA/02_HLA_Histogram_histogram.plotly",
#         details="analysis/cohort_report/HLA/02_details.yaml",
#     params:
#         caption="""caption: 'Histogram of shared HLA-alleles'""",
#         plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0}}),
#     message:
#         "REPORT: creating HLA table for HLA section"
#     group: "cohort_report"
#     shell:
#         """echo "{params.caption}" >> {output.details} && 
#         echo "{params.plot_options}" >> {output.details} &&
#         cidc_wes/modules/scripts/cohort_report/cr_hla_hlaHistogram.py -f {input} -o {output.csv}"""

rule cohort_report_HLA_oncoplot:
    """Generate the HLA histogram plot for the report"""
    input:
        csv="analysis/cohort_report/HLA/02_HLA_Table.dt",
    output:
        csv="analysis/cohort_report/HLA/01_HLA_Oncoplot.stub",
        details="analysis/cohort_report/HLA/01_details.yaml",
    params:
        caption="""caption: 'Oncoplot of shared HLA-alleles'""",
    message:
        "REPORT: creating HLA oncoplot for HLA section"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        touch {output.csv}"""
        #cidc_wes/modules/scripts/cohort_report/cr_hla_hlaOncoplot.py -f {input} -o {output.csv}"""

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
#MISC section
rule cohort_report_tcellextrect:
    """Generate the tcell estimate plot for the report"""
    output:
        stub="analysis/cohort_report/MISC/01_tcellextrect.stub",
        details="analysis/cohort_report/MISC/01_details.yaml",
    params:
        caption="""caption: 'TCell Fraction estimated by TcellExTRECT'""",
    message:
        "REPORT: creating tcell estimate report"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        touch {output.stub}"""
    
rule cohort_report_msisensor2:
    """Generate the msisensor2 plot for the report"""
    output:
        stub="analysis/cohort_report/MISC/02_msisensor2.stub",
        details="analysis/cohort_report/MISC/02_details.yaml",
    params:
        caption="""caption: 'Microsatellite instability estimates'""",
    message:
        "REPORT: creating msisensor report"
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} && 
        touch {output.stub}"""
    
###############################################################################
rule cohort_report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        cohort_report_targets
    params:
        jinja2_template="cidc_wes/report/index.cohort.html",
        report_path = "analysis/cohort_report",
        sections_list=",".join(['data_quality','copy_number_variation', 'somatic_variants', 'HLA', 'MISC']),
        title="WES Cohort Report",
        #sections_list=",".join(['somatic_variants'])
    output:
        "analysis/cohort_report/report.html"
    message:
        "REPORT: Generating WES report"
    group: "cohort_report"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -j {params.jinja2_template} -t "{params.title}" -o {output} && cp -r cidc_wes/report/static {params.report_path}"""

