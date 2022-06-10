#MODULE: wes report module
from yaml import dump as yaml_dump

def report_targets_sansHTML(wildcards):
    """Generates the targets for this module"""
    ls = []
    #Take first element in runs
    run = list(config['runs'].keys())[0]
    #META
    ls.append("analysis/report/config.yaml")
    ls.append("analysis/report/metasheet.csv")
    ls.append("analysis/report/WES_Meta/02_WES_Run_Version.tsv")
    ls.append("analysis/report/WES_Meta/01_WES_Software_Versions.tsv")
    #Data Quality
    ls.append("analysis/report/data_quality/01_mapping_stats.tsv")
    ls.append("analysis/report/data_quality/02_QC_Plots.tsv")
    ls.append("analysis/report/data_quality/03_coverage_statistics.tsv")

    #SOMATIC
    ls.append("analysis/report/somatic_variants/01_somatic_variants_summary.png")

    ls.append("analysis/report/somatic_variants/02_summary_table.csv")
    ls.append("analysis/report/somatic_variants/03_functional_annotation.csv")
    ls.append("analysis/report/somatic_variants/04_SNV_Statistics.csv")
    if not config.get('tumor_only'): #Only run when we have normals
        ls.append("analysis/report/somatic_variants/05_tumor_germline_overlap.tsv")
    ls.append("analysis/report/somatic_variants/06_VAF.png")
    ls.append("analysis/report/somatic_variants/07_lego_plot.png")
    ls.append("analysis/report/somatic_variants/08_lollipop_plots.csv")

    ls.append("analysis/report/somatic_variants/plots/20_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/plots/21_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/plots/22_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/plots/23_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/plots/24_lollipop_plot.png")


    #COPYNUMBER
    if 'copynumber' not in config['skipped_modules'] and not config.get('tumor_only'):
        ls.append("analysis/report/copy_number_variation/01_copy_number_variation_plot.png")
    if 'clonality' not in config['skipped_modules']:
        ls.append("analysis/report/copy_number_variation/02_tumor_clonality.tsv")
    if 'purity' not in config['skipped_modules']:
        ls.append("analysis/report/copy_number_variation/03_tumor_purity.tsv")

    #NEOANTIGEN
    ls.append("analysis/report/neoantigens/01_HLA_Results.tsv")
    if 'neoantigen' not in config['skipped_modules']:
        ls.append("analysis/report/neoantigens/02_neoantigen_list.dt")
    if 'tcellextrect' not in config['skipped_modules']:
        ls.append("analysis/report/neoantigens/03_tcellextrect.csv")
    if 'msisensor2' not in config['skipped_modules']:
        ls.append("analysis/report/neoantigens/04_msisensor2.csv")

    #JSON
    ls.append("analysis/report/json/%s.wes.json" % run)

    return ls

def report_targets(wildcards):
    ls = report_targets_sansHTML(wildcards)
    #REPORT
    ls.append("analysis/report/report.html")
    ls.append("analysis/report.tar.gz")
    return ls

rule report_all:
    input:
        report_targets
    benchmark: "benchmarks/report/report_all.txt"

###############################################################################
#META information
rule report_meta_version:
    """Gather's all the information required for the meta information and
    outputs a csv file
    """
    input:
        config="config.yaml",
        wes_versions="cidc_wes/static/wes_versions.yaml"
    output:
         "analysis/report/WES_Meta/02_WES_Run_Version.tsv"
    message:
        "REPORT: creating WES version table"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report_meta.py -c {input.config} -v {input.wes_versions} -o {output}"""

rule report_meta_software:
    """Gather's all the information required for the meta information and
    outputs a csv file
    """
    input:
        config="config.yaml",
        wes_versions="cidc_wes/static/wes_versions.yaml"
    output:
         "analysis/report/WES_Meta/01_WES_Software_Versions.tsv"
    message:
        "REPORT: creating WES software table"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report_software.py -c {input.config} -v {input.wes_versions} -o {output}"""

###############################################################################

rule report_data_quality_table:
    """Generate the mapping stats table for the report"""
    input:
        "analysis/align/mapping.csv"
    output:
        tsv="analysis/report/data_quality/01_mapping_stats.tsv",
        details="analysis/report/data_quality/01_details.yaml",
    params:
        caption="""caption: 'This table shows the total number reads in each sample, how many of those reads were mapped, and how many are de-duplicated reads.'"""
    message:
        "REPORT: creating mapping stats for data_quality section"
    group: "report"
    shell:
        """echo "{params.caption}" >> {output.details} && cidc_wes/modules/scripts/report_dataQual_mappingStats.py -f {input} -o {output.tsv}"""

#------------------------------------------------------------------------------
def report_data_quality_plotsInputFn(wildcards):
    """Given a run, returns a list of the various sub plots to generate
    """
    ret = []
    run_name = list(config['runs'].keys())[0]
    run = config['runs'][run_name]
    if config.get('tumor_only'): run= run[1:] #if tumor_only drop normal sample
    for sample in run:
        ret.append("analysis/report/data_quality/plots/%s_gcBias.png" % (sample))
        ret.append("analysis/report/data_quality/plots/%s_qualityScore.png" % (sample))
        ret.append("analysis/report/data_quality/plots/%s_qualityByCycle.png" % (sample))
        ret.append("analysis/report/data_quality/plots/%s_insertSize.png" % (sample))
    return ret

rule report_data_quality_gcPlot:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/data_quality/plots/{sample}_gcBias.png"
    params:
        page = 1
    message:
        "REPORT: generating data_quality gc bias plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_data_quality_qualityScore:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/data_quality/plots/{sample}_qualityScore.png"
    params:
        page = 2
    message:
        "REPORT: generating data_quality quality score plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_data_quality_qualityByCycle:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/data_quality/plots/{sample}_qualityByCycle.png"
    params:
        page = 3
    message:
        "REPORT: generating data_quality quality by cycle plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_data_quality_insertSize:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/data_quality/plots/{sample}_insertSize.png"
    params:
        page = 4
    message:
        "REPORT: generating data_quality insert size plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

def report_getTumorNormal(index):
    run = list(config['runs'].keys())[0]
    return config['runs'][run][index]

rule report_data_quality_plots_table:
    """Generate the fastqc plots table for the report
    The trick here is that the plot table shows BOTH samples, tumor and
    normal, for the run--so we need 1. an inpput fn to look up the samples
    and 2. params to generate the correct plots
    """
    input:
        report_data_quality_plotsInputFn
    output:
        tsv="analysis/report/data_quality/02_QC_Plots.tsv",
        details = "analysis/report/data_quality/02_details.yaml",
    params:
        normal= lambda wildcards: "-n %s" % report_getTumorNormal(0) if not config.get('tumor_only') else "",
        tumor = lambda wildcards: report_getTumorNormal(1),
        image_paths = lambda wildcards: "analysis/report/data_quality/plots/",
        sub_caption = """subcaption: 'NOTE: (T) denotes tumor sample; (N) denotes normal sample. a) GC Plot shows the distribution of %GC bases within a 100bp window.  In human, the mean GC content is approx. 40%. b) Quality Score shows the distribution of phred scores. c) Quality by Cycle shows the phred score across the sequencing cycles. d) Insert size shows the distribution of fragment lengths.' """,
    message:
        "REPORT: creating QC plots for data_quality section"
    group: "report"
    shell:
        """echo "{params.sub_caption}" >> {output.details} && cidc_wes/modules/scripts/report_dataQual_plot_table.py {params.normal} -t {params.tumor} -p {params.image_paths} -o {output}"""


rule report_data_quality_coverage:
    """Generate the coverage data table"""
    input:
        "analysis/metrics/all_sample_summaries.txt"
    output:
        tsv="analysis/report/data_quality/03_coverage_statistics.tsv",
        details="analysis/report/data_quality/03_details.yaml"
    params:
        caption="""caption: 'The following table describes the read depth coverage statistics. With the exception of the Total Reads column, which represents the total number of reads in each sample, all numbers represent reads in targeted regions.'"""
    shell:
        """echo "{params.caption}" >> {output.details} && cidc_wes/modules/scripts/report_dataQual_coverage.py -f {input} -o {output.tsv}"""

###############################################################################
def report_maftoolsPlotsInput(wildcards):
    run = list(config['runs'].keys())[0]
    caller = config.get('somatic_caller', "tnscope")
    return "analysis/somatic/%s/%s_%s.output.twist.maf" % (run, run, caller)

rule report_somatic_variants_maftoolsPlots:
    input:
        report_maftoolsPlotsInput
    output:
        summary="analysis/report/somatic_variants/01_somatic_variants_summary.png",
        vaf="analysis/report/somatic_variants/06_VAF.png",
        lolli1="analysis/report/somatic_variants/plots/20_lollipop_plot.png",
        lolli2="analysis/report/somatic_variants/plots/21_lollipop_plot.png",
        lolli3="analysis/report/somatic_variants/plots/22_lollipop_plot.png",
        lolli4="analysis/report/somatic_variants/plots/23_lollipop_plot.png",
        lolli5="analysis/report/somatic_variants/plots/24_lollipop_plot.png",
    params:
        #cap3 = """caption: 'This table summarizes the number of transitions and transversions occuring within the set of SNPs.'"""
    message:
        "REPORT: creating mafTools plot for somatic section"
    group: "report"
    shell:
        """Rscript cidc_wes/modules/scripts/report_somatic_mafPlots.R {input} {output.summary} {output.vaf} {output.lolli1} {output.lolli2} {output.lolli3} {output.lolli4} {output.lolli5}"""

rule report_somatic_variants_lollipop_table:
    """Generate the table of lollipop plots"""
    input:
        lolli01= "analysis/report/somatic_variants/plots/20_lollipop_plot.png",
        lolli02= "analysis/report/somatic_variants/plots/21_lollipop_plot.png",
        lolli03= "analysis/report/somatic_variants/plots/22_lollipop_plot.png",
        lolli04= "analysis/report/somatic_variants/plots/23_lollipop_plot.png",
        lolli05= "analysis/report/somatic_variants/plots/24_lollipop_plot.png",
    output:
        csv = "analysis/report/somatic_variants/08_lollipop_plots.csv",
        details="analysis/report/somatic_variants/08_details.yaml",
    params:
        caption="""caption: 'This table shows lollipop plots for the top 5 cancer driving genes.'""",
        report_path = "analysis/report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cidc_wes/modules/scripts/cohort_report/cr_somatic_lollipop_table.py -f {input.lolli01} -f {input.lolli02} -f {input.lolli03} -f {input.lolli04} -f {input.lolli05} -r {params.report_path} -o {output}"""


def report_somatic_variants_summary_tblInputFn(wildcards):
    ls = []
    caller = config['somatic_caller']
    run = list(config['runs'].keys())[0]
    ls.append("analysis/somatic/%s/%s_%s.mutation_summaries.csv" % (run,run,caller))
    ls.append("analysis/somatic/%s/%s_%s.functional_annot_summaries.csv" % (run,run,caller))
    ls.append("analysis/somatic/%s/%s_%s_somatic_SNV_summaries.csv" % (run, run, caller))
    return ls

rule report_somatic_variants_summary_tbls:
    input:
        report_somatic_variants_summary_tblInputFn
    params:
        cap3 = """caption: 'This table summarizes the number of transitions and transversions occuring within the set of SNPs.'"""
    output:
        csv1 = "analysis/report/somatic_variants/02_summary_table.csv",
        csv2 = "analysis/report/somatic_variants/03_functional_annotation.csv",
        csv3 = "analysis/report/somatic_variants/04_SNV_Statistics.csv",
        details3 = "analysis/report/somatic_variants/04_details.yaml",
    shell:
        """echo "{params.cap3}" >> {output.details3} && cp {input[0]} {output.csv1} && cp {input[1]} {output.csv2} && cp {input[2]} {output.csv3}"""

def report_legoPlotInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    caller = config['somatic_caller']
    return "analysis/somatic/%s/%s_%s.twist.pdf" % (run, run, caller)

rule report_somatic_variants_legoPlot:
    """Add the lego plot to the somatic_variants section of the report"""
    input:
        report_legoPlotInputFn
    output:
        png = "analysis/report/somatic_variants/07_lego_plot.png",
    params:
        page = 1,
    message:
        "REPORT: generating somatic_variants lego plot"
    group: "report"
    shell:
        """Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output.png} {params.page}"""

def report_somatic_variants_germlineCompareInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/germline/%s/%s_vcfcompare.txt" % (run, run)

rule report_somatic_variants_germlineCompare:
    """report germline overlap"""
    input:
        report_somatic_variants_germlineCompareInputFn
    params:
        run = list(config['runs'].keys())[0],
        #targets = center_targets[config.get('cimac_center', 'broad')],
        cap = """caption: 'The table reports the total number of UNFILTERED variants in the tumor sample and normal sample, the number of mutations that they have in common, and their percent overlap.'""",
        sub = """subcaption: 'NOTE: the % overlap was calculated using the number of tumor variants as the denominator'""",
    output:
        tsv = "analysis/report/somatic_variants/05_tumor_germline_overlap.tsv",
        details = "analysis/report/somatic_variants/05_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && echo "{params.sub}" >> {output.details} && cidc_wes/modules/scripts/report_somatic_overlap.py -f {input} -r {params.run} -o {output.tsv}"""

###############################################################################
def report_copynumberPlotInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/clonality/%s/%s_genome_view.pdf" % (run,run)

rule report_copynumberPlot:
    """report tumor copynumberPlot"""
    input:
        report_copynumberPlotInputFn
    params:
        pg = 3,
        subcap = """subcaption: 'Genome-whide visualization of the allele-specific and absolute copy number results, and raw profile of the depth ratio and allele frequency. (ref: https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.html#plots-and-results)'"""
    output:
        png="analysis/report/copy_number_variation/01_copy_number_variation_plot.png",
        details="analysis/report/copy_number_variation/01_details.yaml",
    shell:
        """echo "{params.subcap}" >> {output.details} && Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output.png} {params.pg}"""

def report_copy_number_purityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run)

rule report_copy_number_purity:
    """report tumor purity"""
    input:
        report_copy_number_purityInputFn
    params:
        run = list(config['runs'].keys())[0],
        cap = """caption: 'This table reports the estimated tumor purity and ploidy of the sample.'"""
    output:
        tsv="analysis/report/copy_number_variation/03_tumor_purity.tsv",
        details="analysis/report/copy_number_variation/03_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/report_cnv_purity.py -f {input} -r {params.run} -o {output.tsv}"""

def report_copy_number_clonalityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/clonality/%s/%s_pyclone6.results.summary.tsv" % (run,run)

rule report_copy_number_clonality:
    """report tumor clonality"""
    input:
        report_copy_number_clonalityInputFn
    params:
        run = list(config['runs'].keys())[0],
        cap = """caption: 'This table reports the estimated cancer cell fraction of each cluster.  NOTE: these estimates were based on this single sample's Copy Number Information; for more reliable estimates, multiple samples must be used.'"""
    output:
        tsv="analysis/report/copy_number_variation/02_tumor_clonality.tsv",
        details="analysis/report/copy_number_variation/02_details.yaml",
    shell:
        #"""echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/report_cnv_clonality.py -f {input} -r {params.run} -o {output.tsv}"""
        """echo "{params.cap}" >> {output.details} && cp {input} {output.tsv}"""

###############################################################################
def report_neoantigens_HLAInputFn(wildcards):
    runName = list(config['runs'].keys())[0]
    run = config['runs'][runName]
    normal = run[0]
    tumor = run[1]

    tmp = {}
    if config.get('neoantigen_run_classII'):
        if not config.get('tumor_only'): #Only run when we have normals
            tmp['normal_opti'] = "analysis/optitype/%s/%s_result.tsv" % (normal, normal)
            tmp['normal_hlahd'] = "analysis/hlahd/%s/result/%s_final.result.txt" % (normal, normal)
        tmp['tumor_opti'] = "analysis/optitype/%s/%s_result.tsv" % (tumor, tumor)
        tmp['tumor_hlahd'] = "analysis/hlahd/%s/result/%s_final.result.txt" % (tumor, tumor)
    else:
        #optitype only
        if not config.get('tumor_only'): #Only run when we have normals
            tmp['normal_opti'] = "analysis/optitype/%s/%s_result.tsv" % (normal, normal)
	#tmp['tumor_opti'] = "analysis/optitype/%s/%s_result.tsv" % (tumor, normal)
        tmp['tumor_opti'] = "analysis/optitype/%s/%s_result.tsv" % (tumor, tumor)

    return tmp

def report_getSampleNames():
    """returns the sample names associated with the run.  if tumor_only,
    returns only the tumor sample name"""
    if config.get('tumor_only'):
        return [report_getTumorNormal(1)]
    else:
        return [report_getTumorNormal(0), report_getTumorNormal(1)]

rule report_neoantigens_HLA:
    """report HLA type"""
    input:
        unpack(report_neoantigens_HLAInputFn)
    params:
        #names = ",".join(config['runs'][list(config['runs'].keys())[0]]),
        names = ",".join(report_getSampleNames()),
        cap = """caption: 'This table shows the HLA alleles for both tumor and normal samples.'""",
        #THIS LINE SHOULD BE REPLACED WITH A FUNCTION CALL TO DETERMINE IF XHLA SHOULD BE RUN
        in_files = lambda wildcards,input: "-t %s -u %s" % (input.tumor_opti, input.tumor_hlahd) if config.get('tumor_only', False) else "-n %s -m %s -t %s -u %s" % (input.normal_opti, input.normal_hlahd, input.tumor_opti, input.tumor_hlahd)
    output:
        tsv="analysis/report/neoantigens/01_HLA_Results.tsv",
        details="analysis/report/neoantigens/01_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/report_neoantigens_hla.py {params.in_files} -s {params.names} -o {output.tsv}"""

def report_neoantigens_neoantigen_listInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/neoantigen/%s/%s_neoantigen_table.tsv" % (run, run)

rule report_neoantigens_neoantigen_list:
    """report HLA type"""
    input:
        report_neoantigens_neoantigen_listInputFn
    params:
        cap = """caption: 'This table shows the list of predicted neoantigens.'""",
    output:
        tsv= "analysis/report/neoantigens/02_neoantigen_list.dt",
        details= "analysis/report/neoantigens/02_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cp {input} {output.tsv}"""

def report_neoantigens_tcellextrectInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/tcellextrect/%s/%s_tcellextrect.txt" % (run, run)

rule report_neoantigens_tcellextrect:
    """report T-cell fraction"""
    input:
        report_neoantigens_tcellextrectInputFn
    params:
        #run = lambda wildcards: wildcards.run,
        cap = """caption: 'This table shows the estimated T-cell fraction and associated QC value.'"""
    output:
        table="analysis/report/neoantigens/03_tcellextrect.csv",
	details= "analysis/report/neoantigens/03_details.yaml",
    group: "report"
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/tcellextrect_trimTable.py -f {input} -o {output.table}"""


def report_neoantigens_msisensor2InputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/msisensor2/%s/%s_msisensor2.txt" % (run, run)

rule report_neoantigens_msisensor2:
    """report microsatellite instability"""
    input:
        report_neoantigens_msisensor2InputFn
    params:
        #run = lambda wildcards: wildcards.run,
        cap = """caption: 'This table shows the estimated number and percentage of somatic microsatellite sites.'"""
    output:
        table="analysis/report/neoantigens/04_msisensor2.csv",
        details= "analysis/report/neoantigens/04_details.yaml",
    group: "report"
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/msisensor2_formatTable.py -f {input} -o {output.table}"""


###############################################################################

rule report_json_hla:
    """encode sample hla alleles as json"""
    input:
        unpack(report_neoantigens_HLAInputFn)
    output:
        "analysis/report/json/hla/{run}.hla.json"
    params:
        run = lambda wildcards: wildcards.run,
        in_files = lambda wildcards,input: "-t %s -u %s" % (input.tumor_opti, input.tumor_hlahd) if config.get('tumor_only', False) else "-n %s -m %s -t %s -u %s" % (input.normal_opti, input.normal_hlahd, input.tumor_opti, input.tumor_hlahd)
    group: "report"
    shell:
        "cidc_wes/modules/scripts/json_report_hla.py -r {params.run} {params.in_files} -o {output}"

###############################################################################
rule report_copy_runInfoFiles:
    input:
        config="config.yaml",
        metasheet= "metasheet.csv",
    output:
        conf="analysis/report/config.yaml",
        meta="analysis/report/metasheet.csv",
    shell:
        """cp {input.config} {output.conf} &&
        cp {input.metasheet} {output.meta}"""

def getJsonFiles(wildcards):
    """Will return a list of the run's json files"""
    run = wildcards.run
    caller = config.get("somatic_caller", "tnscope")
    tmp = {}
    tmp['mapping'] = "analysis/report/json/align/%s.mapping.json" % run
    tmp['coverage'] = "analysis/report/json/coverage/%s.coverage.json" % run
    tmp['gc_content'] = "analysis/report/json/gc_content/%s.gc_content.json" % run
    tmp['insert_size'] = "analysis/report/json/insert_size/%s.insert_size.json" % run
    tmp['mean_quality'] = "analysis/report/json/mean_quality/%s.mean_quality.json" % run
    tmp['hla'] = "analysis/report/json/hla/%s.hla.json" % run
    tmp['somatic'] = "analysis/report/json/somatic/%s_%s.somatic.json" % (run, caller)
    tmp['neoantigen'] = "analysis/report/json/neoantigen/%s.neoantigen.json" % run
    tmp['tcellextrect'] = "analysis/report/json/tcellextrect/%s.tcellextrect.json" % run
    tmp['msisensor2'] = "analysis/report/json/msisensor2/%s.msisensor2.json" % run
    tmp['purity'] = "analysis/report/json/purity/%s.purity.json" % run
    tmp['clonality'] = "analysis/report/json/clonality/%s.clonality.json" % run
    tmp['copynumber'] = "analysis/report/json/copynumber/%s.copynumber.json" % run

    # remove optional modules that are not being run
    for module in ['purity', 'clonality', 'neoantigen', 'msisensor2', 'tcellextrect']:
        if module in config['skipped_modules']:
            tmp.pop(module)

    return tmp

#LEN: REVISE!--explicit param calls!
def buildJsonParams(file_dict):
    '''maps input jsons to thier commandline arguments and returns a concatenated string'''
    run = list(config['runs'].keys())[0]
    caller = config.get("somatic_caller", "tnscope")
    arg_dict = {'mapping': '-m', 'coverage': '-c', 'gc_content': '-g',
                'insert_size': '-i', 'mean_quality':'-q', 'hla': '-j',
                'somatic': '-s', 'clonality': '-l', 'purity': '-p', 'neoantigen': '-n',
                'tcellextrect': '-t', 'msisensor2': '-e', 'copynumber': '-v',
    }

    ret = ''
    for module in file_dict.keys():
        ret = ret + arg_dict[module] + ' ' + file_dict[module] + ' '

    return ret


#LEN: REVISE!--explicit param calls!
rule report_generate_json:
    input:
        unpack(getJsonFiles)
    params:
        run = lambda wildcards: wildcards.run,
        in_files = lambda wildcards,input: buildJsonParams(input)
    output:
        #NOTE: CANNOT name this {run}.json otherwise snakemake will have
        #trouble ressolving the wildcard
        "analysis/report/json/{run}.wes.json"
    shell:
        """cidc_wes/modules/scripts/json_stitcher.py -r {params.run} {params.in_files} -o {output}"""


###############################################################################

def getSections():
    '''creates string of tabs needed for the report based on which modules are run'''
    sections = ['WES_Meta','data_quality', 'copy_number_variation','somatic_variants','neoantigens']
    not_copynumber_modules = ('copynumber' in config['skipped_modules']) and ('clonality' in config['skipped_modules']) and ('purity' in config['skipped_modules'])

    if not_copynumber_modules or config.get('tumor_only'):
        sections.remove('copy_number_variation')

    sections_str = ','.join(sections)
    return sections_str
    


rule report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report_targets_sansHTML
    params:
        jinja2_template="cidc_wes/report/index.sample.html",
        report_path = "analysis/report",
	sections_list = getSections()
    output:
        "analysis/report/report.html"
    message:
        "REPORT: Generating WES report"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report.py -d {params.report_path} -s {params.sections_list} -j {params.jinja2_template} -o {output} && cp -r cidc_wes/report/static {params.report_path}"""

rule report_gzipReport:
    input:
        "analysis/report/report.html"
    output:
        "analysis/report.tar.gz"
    message: "REPORT: Zipping up report directory"
    group: "report"
    shell:
        "tar -c analysis/report | gzip > {output}"
