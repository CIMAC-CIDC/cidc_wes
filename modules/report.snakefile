#MODULE: wes report module 
from yaml import dump as yaml_dump

def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    #META
    ls.append("analysis/report/config.yaml")
    ls.append("analysis/report/metasheet.csv")
    ls.append("analysis/report/wes_meta/02_wes_run_version.tsv")
    ls.append("analysis/report/wes_meta/01_wes_software_versions.tsv")
    #Data Quality
    ls.append("analysis/report/data_quality/01_mapping_stats.tsv")
    ls.append("analysis/report/data_quality/02_qc_plots.tsv")
    ls.append("analysis/report/data_quality/03_coverage_statistics.tsv")
    #SOMATIC
    ls.append("analysis/report/somatic_variants/03_summary_table.csv")
    ls.append("analysis/report/somatic_variants/04_functional_annotation.csv")
    ls.append("analysis/report/somatic_variants/05_SNV_statistics.csv")

    #somatic - MAFtools plots---REALLY need to reorganize the somatic section
    ls.append("analysis/report/somatic_variants/01_somatic_variants_summary.png")
    ls.append("analysis/report/somatic_variants/02_vaf.png")

    ls.append("analysis/report/somatic_variants/09_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/10_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/11_lollipop_plot.png")
    ls.append("analysis/report/somatic_variants/12_lollipop_plot.png")

    ls.append("analysis/report/somatic_variants/07_tumor_mutational_burden.tsv")
    #COPYNUMBER
    ls.append("analysis/report/copy_number/01_copynumber_plot.png")
    ls.append("analysis/report/copy_number/02_tumor_clonality.tsv")
    ls.append("analysis/report/copy_number/03_tumor_purity.tsv")
    
    #NEOANTIGEN
    ls.append("analysis/report/neoantigens/01_HLA_results.tsv")
    ls.append("analysis/report/neoantigens/02_neoantigen_list.tsv")
        
    for run in config['runs']:
        ls.append("analysis/report/somatic_variants/08_%s_lego_plot.png" % run)
    return ls

rule report_all:
    input:
        "analysis/report/report.html",
        "analysis/report.tar.gz"

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
         "analysis/report/wes_meta/02_wes_run_version.tsv"
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
         "analysis/report/wes_meta/01_wes_software_versions.tsv"
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
        tsv="analysis/report/data_quality/02_qc_plots.tsv",
        details = "analysis/report/data_quality/02_details.yaml",
    params:
        normal= lambda wildcards: report_getTumorNormal(0),
        tumor = lambda wildcards: report_getTumorNormal(1),
        image_paths = lambda wildcards: "analysis/report/data_quality/plots/",
        sub_caption = """subcaption: 'NOTE: (T) denotes tumor sample; (N) denotes normal sample. a) GC Plot shows the distribution of %GC bases within a 100bp window.  In human, the mean GC content is approx. 40%. b) Quality Score shows the distribution of phred scores. c) Quality by Cycle shows the phred score across the sequencing cycles. d) Insert size shows the distribution of fragment lengths.' """,
    message:
        "REPORT: creating QC plots for data_quality section"
    group: "report"
    shell:
        """echo "{params.sub_caption}" >> {output.details} && cidc_wes/modules/scripts/report_dataQual_plot_table.py -n {params.normal} -t {params.tumor} -p {params.image_paths} -o {output}"""

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
    return "analysis/somatic/%s/%s_%s.filter.maf" % (run, run, caller)

rule report_somatic_variants_maftoolsPlots:
    input:
        report_maftoolsPlotsInput
    output:
        summary="analysis/report/somatic_variants/01_somatic_variants_summary.png",
        vaf="analysis/report/somatic_variants/02_vaf.png",
        lolli1="analysis/report/somatic_variants/09_lollipop_plot.png",
        lolli2="analysis/report/somatic_variants/10_lollipop_plot.png",
        lolli3="analysis/report/somatic_variants/11_lollipop_plot.png",
        lolli4="analysis/report/somatic_variants/12_lollipop_plot.png",
    params:
        #cap3 = """caption: 'This table summarizes the number of transitions and transversions occuring within the set of SNPs.'"""
    message:
        "REPORT: creating mafTools plot for somatic section"
    group: "report"
    shell:
        """Rscript cidc_wes/modules/scripts/report_somatic_mafPlots.R {input} {output.summary} {output.vaf} {output.lolli1} {output.lolli2} {output.lolli3} {output.lolli4}"""


def report_somatic_variants_summary_tblInputFn(wildcards):
    ls = []
    caller = config['somatic_caller']
    run = list(config['runs'].keys())[0]
    ls.append("analysis/somatic/somatic_mutation_summaries.%s.csv" % caller)
    ls.append("analysis/somatic/somatic_functional_annot_summaries.%s.csv" % caller)
    ls.append("analysis/somatic/%s/%s_%s_somatic_SNV_summaries.csv" % (run, run, caller))
    return ls

rule report_somatic_variants_summary_tbls:
    input:
        report_somatic_variants_summary_tblInputFn
    params:
        cap3 = """caption: 'This table summarizes the number of transitions and transversions occuring within the set of SNPs.'"""
    output:
        csv1 = "analysis/report/somatic_variants/03_summary_table.csv",
        csv2 = "analysis/report/somatic_variants/04_functional_annotation.csv",
        csv3 = "analysis/report/somatic_variants/05_SNV_statistics.csv",
        details3 = "analysis/report/somatic_variants/04_details.yaml",
    shell:
        """echo "{params.cap3}" >> {output.details3} && cp {input[0]} {output.csv1} && cp {input[1]} {output.csv2} && cp {input[2]} {output.csv3}"""

def report_legoPlotInputFn(wildcards):
    run = wildcards.run
    caller = config['somatic_caller']
    return "analysis/somatic/%s/%s_%s.filter.pdf" % (run, run, caller)

rule report_somatic_variants_legoPlot:
    """Add the lego plot to the somatic_variants section of the report"""
    input:
        report_legoPlotInputFn
    output:
        png = "analysis/report/somatic_variants/08_{run}_lego_plot.png",
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
        cap = """caption: 'This table reports the tumor mutational burden (TMB) as well as the total mutational load of the normal sample, the number of mutations that they have in common, and their percent overlap.'""",
        sub = """subcaption: 'NOTE: the % overlap was calculated using the Tumor TMB as the denominator'""",
    output:
        tsv = "analysis/report/somatic_variants/07_tumor_mutational_burden.tsv",
        details = "analysis/report/somatic_variants/07_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && echo "{params.sub}" >> {output.details} && cidc_wes/modules/scripts/report_somatic_tmb.py -f {input} -r {params.run} -o {output.tsv}"""

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
        png="analysis/report/copy_number/01_copynumber_plot.png",
        details="analysis/report/copy_number/01_details.yaml",
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
        cap = """caption: 'This table reports the estimated tumor purity, ploidy, and diploid log ratio of the sample.'"""
    output:
        tsv="analysis/report/copy_number/03_tumor_purity.tsv",
        details="analysis/report/copy_number/03_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/report_cnv_purity.py -f {input} -r {params.run} -o {output.tsv}"""

def report_copy_number_clonalityInputFn(wildcards):
    run = list(config['runs'].keys())[0]
    return "analysis/clonality/%s/%s_table.tsv" % (run,run)

rule report_copy_number_clonality:
    """report tumor clonality"""
    input:
        report_copy_number_clonalityInputFn
    params:
        run = list(config['runs'].keys())[0],
        cap = """caption: 'This table reports the estimated tumor clonaltiy of the sample.'"""
    output:
        tsv="analysis/report/copy_number/02_tumor_clonality.tsv",
        details="analysis/report/copy_number/02_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/report_cnv_clonality.py -f {input} -r {params.run} -o {output.tsv}"""

###############################################################################
def report_neoantigens_HLAInputFn(wildcards):
    runName = list(config['runs'].keys())[0]
    run = config['runs'][runName]
    normal = run[0]
    tumor = run[1]

    tmp = {}
    if config.get('neoantigen_run_classII'):
        tmp['normal_opti'] = "analysis/optitype/%s/%s_result.tsv" % (normal, normal)
        tmp['normal_xhla'] = "analysis/xhla/%s/report-%s-hla.json" % (normal, normal)
        tmp['tumor_opti'] = "analysis/optitype/%s/%s_result.tsv" % (tumor, tumor)
        tmp['tumor_xhla'] = "analysis/xhla/%s/report-%s-hla.json" % (tumor, tumor)
    else:
        #optitype only
        tmp['normal_opti'] = "analysis/optitype/%s/%s_result.tsv" % (normal, normal)
        tmp['tumor_opti'] = "analysis/optitype/%s/%s_result.tsv" % (tumor, normal)
    return tmp

rule report_neoantigens_HLA:
    """report HLA type"""
    input:
        unpack(report_neoantigens_HLAInputFn)
    params:
        names = ",".join(config['runs'][list(config['runs'].keys())[0]]),
        cap = """caption: 'This table shows the HLA alleles for both tumor and normal samples.'""",
    output:
        tsv="analysis/report/neoantigens/01_HLA_results.tsv",
        details="analysis/report/neoantigens/01_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cidc_wes/modules/scripts/report_neoantigens_hla.py -n {input.normal_opti} -m {input.normal_xhla} -t {input.tumor_opti} -u {input.tumor_xhla} -s {params.names} -o {output.tsv}"""

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
        tsv= "analysis/report/neoantigens/02_neoantigen_list.tsv",
        details= "analysis/report/neoantigens/02_details.yaml",
    shell:
        """echo "{params.cap}" >> {output.details} && cp {input} {output.tsv}"""
###############################################################################

rule report_json_hla:
    """encode sample hla alleles as json"""
    input:
        unpack(report_neoantigens_HLAInputFn)
    output:
        "analysis/report/json/hla/{run}.hla.json"
    params:
        run = lambda wildcards: wildcards.run
    group: "report"
    shell:
        "cidc_wes/modules/scripts/json_report_hla.py -r {params.run} -n {input.normal_opti} -m {input.normal_xhla} -t {input.tumor_opti} -u {input.tumor_xhla} -o {output}"

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
    tmp['purity'] = "analysis/report/json/purity/%s.purity.json" % run
    tmp['somatic'] = "analysis/report/json/somatic/%s_%s.filtered_maf.json" % (run, caller)
    tmp['clonality'] = "analysis/report/json/clonality/%s.clonality.json" % run
    return tmp

rule report_generate_json:
    input:
        unpack(getJsonFiles)
    params:
        run = lambda wildcards: wildcards.run
    output:
        #NOTE: CANNOT name this {run}.json otherwise snakemake will have
        #trouble ressolving the wildcard
        "analysis/report/json/{run}.wes.json"
    shell:
        """cidc_wes/modules/scripts/json_stitcher.py -r {params.run} -m {input.mapping} -c {input.coverage} -g {input.gc_content} -i {input.insert_size} -q {input.mean_quality} -j {input.hla} -p {input.purity} -s {input.somatic} -t {input.clonality} -o {output}"""

###############################################################################
rule report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report_targets
    params:
        jinja2_template="cidc_wes/report/index.sample.html",
        report_path = "analysis/report",
        sections_list=",".join(['wes_meta','data_quality', 'copy_number','somatic_variants','neoantigens'])
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
