#MODULE: wes report module 

def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append("analysis/report/wes_meta.html")
    ls.append("analysis/report/wes_level1.html")
    ls.append("analysis/report/wes_level2.html")
    ls.append("analysis/report/wes_level3.html")
    ls.append("analysis/report/static/done.txt")
    #ls.append("analysis/report/wes_images/...")
    for sample in config['samples']:
        ls.append("analysis/report/wes_images/align/%s/%s_gcBias.png" % (sample,sample))
        ls.append("analysis/report/wes_images/align/%s/%s_qualityScore.png" % (sample,sample))
        ls.append("analysis/report/wes_images/align/%s/%s_qualityByCycle.png" % (sample,sample))
        ls.append("analysis/report/wes_images/align/%s/%s_insertSize.png" % (sample,sample))
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

rule report_level3:
    """Generate wes_level3.html"""
    input:
        #NOTE: need to ensure that this runs AFTER everything is generated!
        config="config.yaml"
    output:
         "analysis/report/wes_level3.html"
    message:
        "REPORT: creating wes_level3.html"
    group: "report"
    shell:
        """cidc_wes/modules/scripts/report_level3.py -c {input} -o {output}"""

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
