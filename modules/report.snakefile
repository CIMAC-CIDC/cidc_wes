#MODULE: wes report module 

def pvacseq_plot_inputfn(wildcards):
    """Will return analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.condensed.ranked.addSample.tsv, but will need to derefernce tumor
    USES neoantigen_getTumor fn from neoantigen.snakefile
    """
    ls = []
    run = wildcards.run
    tumor = neoantigen_getTumor(wildcards)[0]
    ls.append("analysis/neoantigen/%s/MHC_Class_I/%s.filtered.condensed.ranked.addSample.tsv" % (run,tumor))
    return ls
    
def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append("analysis/report/wes_meta.html")
    ls.append("analysis/report/wes_level1.html")
    ls.append("analysis/report/wes_level2.html")
    ls.append("analysis/report/static/done.txt")
    ls.append("analysis/report/wes_images/align/mapping.png")
    for sample in config['samples']:
        ls.append("analysis/report/wes_images/align/%s/%s_gcBias.png" % (sample,sample))
        ls.append("analysis/report/wes_images/align/%s/%s_qualityScore.png" % (sample,sample))
        ls.append("analysis/report/wes_images/align/%s/%s_qualityByCycle.png" % (sample,sample))
        ls.append("analysis/report/wes_images/align/%s/%s_insertSize.png" % (sample,sample))
        
    for run in config['runs']:
        ls.append("analysis/report/wes_images/somatic/%s/%s_%s.legoPlot.png" % (run, run, config['somatic_caller']))
        #pvacseq images
        ls.append("analysis/report/wes_images/neoantigen/%s/HLA_epitopes_fraction_plot.png" % run)
        ls.append("analysis/report/wes_images/neoantigen/%s/Patient_count_epitopes_plot.png" % run)
        ls.append("analysis/report/wes_images/neoantigen/%s/epitopes_affinity_plot.png" % run)
        #copynumber
        ls.append("analysis/report/wes_images/copynumber/%s.%s/circos.png" % (run, config['somatic_caller']))

        #TEST if clonality was run or not
        clonality_density = "analysis/clonality/%s/%s_plot.density.pdf" % (run,run)
        if os.path.exists(clonality_density):
            #COPY over the clonality plots
            ls.append("analysis/report/wes_images/clonality/%s/%s_plot.density.png" % (run,run))
            ls.append("analysis/report/wes_images/clonality/%s/%s_plot.scatter.png" % (run,run))
            ls.append("analysis/report/wes_images/clonality/%s/%s_plot.coordinates.png" % (run,run))
    return ls

rule report_all:
    input:
        report_targets

rule report_meta:
    """Generate wes_meta.html"""
    input:
        #NOTE: need to ensure that this runs AFTER everything is generated!
        config="config.yaml"
    output:
         "analysis/report/wes_meta.html"
    message:
        "REPORT: creating wes_meta.html"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report_meta.py -c {input} -o {output}"""

rule report_level1_gcBiasPlot:
    """Generate gcBiasPlot"""
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/wes_images/align/{sample}/{sample}_gcBias.png"
    params:
        page = 1
    message:
        "REPORT: generating wes level1 gc bias plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level1_qualityScore:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/wes_images/align/{sample}/{sample}_qualityScore.png"
    params:
        page = 2
    message:
        "REPORT: generating wes level1 quality score plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level1_qualityByCycle:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/wes_images/align/{sample}/{sample}_qualityByCycle.png"
    params:
        page = 3
    message:
        "REPORT: generating wes level1 quality by cycle plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level1_insertSize:
    input:
        "analysis/metrics/{sample}/{sample}_metrics.pdf"
    output:
        "analysis/report/wes_images/align/{sample}/{sample}_insertSize.png"
    params:
        page = 4
    message:
        "REPORT: generating wes level1 insert size plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level1_somatic_legoPlot:
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.pdf"
    output:
        "analysis/report/wes_images/somatic/{run}/{run}_{caller}.legoPlot.png"
    params:
        page = 1
    message:
        "REPORT: generating wes level1 lego plot"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level1:
    """Generate wes_level1.html"""
    input:
        #NOTE: need to ensure that this runs AFTER everything is generated!
        config="config.yaml"
    output:
         "analysis/report/wes_level1.html"
    message:
        "REPORT: creating wes_level1.html"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report_level1.py -c {input} -o {output}"""

rule report_level2:
    """Generate wes_level2.html"""
    input:
        #NOTE: need to ensure that this runs AFTER everything is generated!
        config="config.yaml"
    output:
         "analysis/report/wes_level2.html"
    message:
        "REPORT: creating wes_level2.html"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report_level2.py -c {input} -o {output}"""

# rule report_level3:
#     """Generate wes_level3.html"""
#     input:
#         #NOTE: need to ensure that this runs AFTER everything is generated!
#         config="config.yaml"
#     output:
#          "analysis/report/wes_level3.html"
#     message:
#         "REPORT: creating wes_level3.html"
#     group: "report"
#     shell:
#         """cidc_wes/modules/scripts/report_level3.py -c {input} -o {output}"""

rule report_cp_static:
    """Copy cidc_wes/reprt/static to analysis/report/static
    HACK: need a 'done' file to indicate the job worked"""
    input:
    output:
         "analysis/report/static/done.txt"
    message:
        "REPORT: copying static files"
    group: "report"
    shell:
        "cp -r cidc_wes/report/static/ analysis/report/ && touch {output}"

rule pvacseq_plot:
    """Plot the three pvacseq images"""
    input:
        pvacseq_plot_inputfn
    output:
        hla="analysis/report/wes_images/neoantigen/{run}/HLA_epitopes_fraction_plot.png",
        patient="analysis/report/wes_images/neoantigen/{run}/Patient_count_epitopes_plot.png",
        epitope="analysis/report/wes_images/neoantigen/{run}/epitopes_affinity_plot.png",
    params:
        outdir = lambda wildcards: "analysis/report/wes_images/neoantigen/%s/" % wildcards.run
    message:
        "REPORT: generating pvacseq images"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/pvacseq_plot.R -i {input} -o {params.outdir}"

rule mapping_plot:
    """Plot the mapping stats"""
    input:
        "analysis/align/mapping.csv"
    output:
        "analysis/report/wes_images/align/mapping.png"
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/map_stats.R {input} {output}"

rule copynumber_circos_plot:
    """Generates the circos plot for a run"""
    input:
        cnv="analysis/copynumber/{run}/{run}_cnvcalls.circos.txt",
        indel="analysis/somatic/{run}/{run}_{caller}.indel.circos.txt",
        snp="analysis/somatic/{run}/{run}_{caller}.snp.circos.txt",
    params:
        output_dir = lambda wildcards, input, output: "/".join(output[0].split("/")[:-1])
    output:
        #NOTE: because the etc script generates a file called circos.png
        #we put the caller in the dir name--not pretty but i don't want to
        #dynamically generatethe etc/circos.conf
        "analysis/report/wes_images/copynumber/{run}.{caller}/circos.png"
    group: "report"
    shell:
        #FOR circos we need to do the following:
        #0. make a data sub-dir in the {output_dir}
        #1. copy/link in the input file as {output_dir}/data/data.*.txt
        #2. copy cidc_wes/static/circos/etc/ into {output_dir}
        #3. run circos (in that directory)
        """mkdir -p {params.output_dir}/data && \
        cp {input.cnv} {params.output_dir}/data/data.cnv.txt && \
        cp {input.indel} {params.output_dir}/data/data.indel.txt && \
        cp {input.snp} {params.output_dir}/data/data.snp.txt && \
        cp -r cidc_wes/static/circos/etc {params.output_dir} && \
        cd {params.output_dir} && circos"""

rule report_level2_density_plot:
    "convert the density.pdf to png"
     input:
         "analysis/clonality/{run}/{run}_plot.density.pdf",
    output:
        "analysis/report/wes_images/clonality/{run}/{run}_plot.density.png",
    params: page=1
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level2_scatter_plot:
    "convert the density.pdf to png"
     input:
         "analysis/clonality/{run}/{run}_plot.scatter.pdf",
    output:
        "analysis/report/wes_images/clonality/{run}/{run}_plot.scatter.png",
    params: page=1
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"

rule report_level2_coordinates_plot:
    "convert the density.pdf to png"
     input:
         "analysis/clonality/{run}/{run}_plot.coordinates.pdf",
    output:
        "analysis/report/wes_images/clonality/{run}/{run}_plot.coordinates.png"
    params: page=1
    group: "report"
    shell:
        "Rscript cidc_wes/modules/scripts/wes_pdf2png.R {input} {output} {params.page}"
