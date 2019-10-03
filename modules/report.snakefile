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
