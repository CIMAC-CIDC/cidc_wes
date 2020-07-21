#MODULE: wes report2 module 
from yaml import dump as yaml_dump

def report2_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #META
    ls.append("analysis/report2/wes_meta/02_wes_run_version.tsv")
    ls.append("analysis/report2/wes_meta/01_wes_software_versions.tsv")
    #Data Quality
    ls.append("analysis/report2/data_quality/01_mapping_stats.tsv")
    ls.append("analysis/report2/data_quality/02_qc_plots.tsv")
    ls.append("analysis/report2/data_quality/03_coverage_statistics.tsv")
    #SOMATIC
    ls.append("analysis/report2/somatic_variants/01_summary_table.csv")
    ls.append("analysis/report2/somatic_variants/02_functional_annotation.csv")
    ls.append("analysis/report2/somatic_variants/03_SNV_statistics.csv")
    ls.append("analysis/report2/somatic_variants/04_tumor_mutational_burden.tsv")
    #COPYNUMBER
    ls.append("analysis/report2/copy_number/01_copynumber_plot.png")
    ls.append("analysis/report2/copy_number/02_tumor_clonality.tsv")
    ls.append("analysis/report2/copy_number/03_tumor_purity.tsv")
    
    #NEOANTIGEN
    ls.append("analysis/report2/neoantigens/01_HLA_results.tsv")
    ls.append("analysis/report2/neoantigens/02_neoantigen_list.tsv")
    
    for run in config['runs']:
        ls.append("analysis/report2/somatic_variants/05_%s_lego_plot.png" % run)
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
         "analysis/report2/wes_meta/02_wes_run_version.tsv"
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
         "analysis/report2/wes_meta/01_wes_software_versions.tsv"
    message:
        "REPORT: creating WES software table"
    group: "report2"
    shell:
        """cidc_wes/modules/scripts/report_software.py -c {input.config} -v {input.wes_versions} -o {output}"""

###############################################################################

rule report2_data_quality_table:
    """Generate the mapping stats table for the report"""
    input:
        "analysis/align/mapping.csv"
    output:
        tsv="analysis/report2/data_quality/01_mapping_stats.tsv",
        cap="analysis/report2/data_quality/01_caption.txt",
    params:
        caption="This table shows the total number reads in each sample, how many of those reads were mapped, and how many are de-duplicated reads."
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
    run_name = list(config['runs'].keys())[0]
    run = config['runs'][run_name]
    for sample in run:
        ret.append("analysis/report2/data_quality/plots/%s_gcBias.png" % (sample))
        ret.append("analysis/report2/data_quality/plots/%s_qualityScore.png" % (sample))
        ret.append("analysis/report2/data_quality/plots/%s_qualityByCycle.png" % (sample))
        ret.append("analysis/report2/data_quality/plots/%s_insertSize.png" % (sample))
    return ret

rule report2_data_quality_gcPlot:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report2/data_quality/plots/{sample}_gcBias.png"
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
        "analysis/report2/data_quality/plots/{sample}_qualityScore.png"
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
        "analysis/report2/data_quality/plots/{sample}_qualityByCycle.png"
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
        "analysis/report2/data_quality/plots/{sample}_insertSize.png"
    params:
        page = 4
    message:
        "REPORT: generating data_quality insert size plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

def report2_getTumorNormal(index):
    run = list(config['runs'].keys())[0]
    return config['runs'][run][index]

rule report2_data_quality_plots_table:
    """Generate the fastqc plots table for the report
    The trick here is that the plot table shows BOTH samples, tumor and 
    normal, for the run--so we need 1. an inpput fn to look up the samples
    and 2. params to generate the correct plots
    """
    input:
        report2_data_quality_plotsInputFn
    output:
        tsv="analysis/report2/data_quality/02_qc_plots.tsv",
        sub_cap = "analysis/report2/data_quality/02_subcaption.txt",
        #cap=...
    params:
        normal= lambda wildcards: report2_getTumorNormal(0),
        tumor = lambda wildcards: report2_getTumorNormal(1),
        image_paths = lambda wildcards: "analysis/report2/data_quality/plots/",
        sub_caption = "NOTE: (T) denotes tumor sample; (N) denotes normal sample. a) GC Plot shows the distribution of %GC bases within a 100bp window.  In human, the mean GC content is approx. 40%. b) Quality Score shows the distribution of phred scores. c) Quality by Cycle shows the phred score across the sequencing cycles. d) Insert size shows the distribution of fragment lengths.",
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
        tsv="analysis/report2/data_quality/03_coverage_statistics.tsv",
        cap="analysis/report2/data_quality/03_caption.txt"
    params:
        caption="""The following table describes the read depth coverage statistics. With the exception of the Total Reads column, which represents the total number of reads in each sample, all numbers represent reads in targeted regions."""
    shell:
        """echo "{params.caption}" > {output.cap} && cidc_wes/modules/scripts/report_dataQual_coverage.py -f {input} -o {output.tsv}"""

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
    params:
        cap3 = "This table summarizes the number of transitions and transversions occuring within the set of SNPs."
    output:
        csv1 = "analysis/report2/somatic_variants/01_summary_table.csv",
        csv2 = "analysis/report2/somatic_variants/02_functional_annotation.csv",
        csv3 = "analysis/report2/somatic_variants/03_SNV_statistics.csv",
        cap3 = "analysis/report2/somatic_variants/03_caption.txt",
    shell:
        """echo "{params.cap3}" > {output.cap3} && cp {input[0]} {output.csv1} && cp {input[1]} {output.csv2} && cp {input[2]} {output.csv3}"""

def report_legoPlotInputFn(wildcards):
    run = wildcards.run
    caller = config['somatic_caller']
    return "analysis/somatic/%s/%s_%s.filter.pdf" % (run, run, caller)

rule report2_somatic_variants_legoPlot:
    """Add the lego plot to the somatic_variants section of the report"""
    input:
        report_legoPlotInputFn
    output:
        png = "analysis/report2/somatic_variants/05_{run}_lego_plot.png",
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
        cap = "This table reports the tumor mutational burden (TMB) as well as the total mutational load of the normal sample, the number of mutations that they have in common, and their percent overlap. ",
        sub = "NOTE: the % overlap was calculated using the Tumor TMB as the denominator",
    output:
        tsv = "analysis/report2/somatic_variants/04_tumor_mutational_burden.tsv",
        cap = "analysis/report2/somatic_variants/04_caption.txt",
        subcap = "analysis/report2/somatic_variants/04_subcaption.txt",
    shell:
        """echo "{params.cap}" > {output.cap} && echo "{params.sub}" > {output.subcap} && cidc_wes/modules/scripts/report_somatic_tmb.py -f {input} -r {params.run} -o {output.tsv}"""

###############################################################################
def report2_copynumberPlotInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/clonality/%s/%s_genome_view.pdf" % (run,run)

rule report2_copynumberPlot:
    """report tumor copynumberPlot"""
    input:
        report2_copynumberPlotInputFn
    params:
        pg = 3,
        subcap = "Genome-whide visualization of the allele-specific and absolute copy number results, and raw profile of the depth ratio and allele frequency. (ref: https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html#plots-and-results)"
    output:
        png="analysis/report2/copy_number/01_copynumber_plot.png",
        subcap="analysis/report2/copy_number/01_subcaption.txt",
    shell:
        """echo "{params.subcap}" > {output.subcap} && Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output.png} {params.pg}"""

def report2_copy_number_purityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run)

rule report2_copy_number_purity:
    """report tumor purity"""
    input:
        report2_copy_number_purityInputFn
    params:
        run = list(config['runs'].keys())[0],
        cap = "This table reports the estimated tumor purity, ploidy, and diploid log ratio of the sample."
    output:
        tsv="analysis/report2/copy_number/03_tumor_purity.tsv",
        cap="analysis/report2/copy_number/03_caption.txt",
    shell:
        """echo "{params.cap}" > {output.cap} && cidc_wes/modules/scripts/report_cnv_purity.py -f {input} -r {params.run} -o {output.tsv}"""

def report2_copy_number_clonalityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/clonality/%s/%s_table.tsv" % (run,run)

rule report2_copy_number_clonality:
    """report tumor clonality"""
    input:
        report2_copy_number_clonalityInputFn
    params:
        run = list(config['runs'].keys())[0],
        cap = "This table reports the estimated tumor clonaltiy of the sample."
    output:
        tsv="analysis/report2/copy_number/02_tumor_clonality.tsv",
        cap="analysis/report2/copy_number/02_caption.txt",
    shell:
        """echo "{params.cap}" > {output.cap} && cidc_wes/modules/scripts/report_cnv_clonality.py -f {input} -r {params.run} -o {output.tsv}"""

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
        cap = "This table shows the HLA alleles for both tumor and normal samples.",
    output:
        tsv="analysis/report2/neoantigens/01_HLA_results.tsv",
        cap="analysis/report2/neoantigens/01_caption.txt",
    shell:
        """echo "{params.cap}" > {output.cap} && cidc_wes/modules/scripts/report_neoantigens_hla.py -n {params.normal} -t {params.tumor} -s {params.names} -o {output.tsv}"""

def report2_neoantigens_neoantigen_listInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/neoantigen/%s/%s_neoantigen_table.tsv" % (run, run)

rule report2_neoantigens_neoantigen_list:
    """report HLA type"""
    input:
        report2_neoantigens_neoantigen_listInputFn
    params:
        cap = "This table shows the list of predicted neoantigens.",
    output:
        tsv= "analysis/report2/neoantigens/02_neoantigen_list.tsv",
        cap= "analysis/report2/neoantigens/02_caption.txt",
    shell:
        """echo "{params.cap}" > {output.cap} && cut -f 1,3,4,5,6,7,8,8,10 {input} > {output.tsv}"""

###############################################################################
rule report2_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report2_targets
    params:
        report_path = "analysis/report2",
        sections_list=",".join(['wes_meta','data_quality', 'copy_number','somatic_variants','neoantigens'])
    output:
        "analysis/report2/report.html"
    message:
        "REPORT: Generating WES report"
    group: "report2"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -o {output} && cp -r cidc_wes/report2/static {params.report_path}"""

