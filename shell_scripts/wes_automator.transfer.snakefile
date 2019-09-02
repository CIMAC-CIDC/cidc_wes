import glob

def all_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append("wes_tx.analysis.txt")
    ls.append("wes_tx.benchmarks.txt")
    ls.append("wes_tx.src.txt")
    ls.append("wes_tx.config_meta.txt")
    ls.append("wes_tx.nohup.txt")
    return ls

configfile: "config.yaml"

rule target:
    input:
        all_targets

rule transfer_analysis:
    input:
        "analysis/align/mapping.csv" #stub file b/c snkmk doesn't allow dirs
    params:
        transfer_bucket=config['transfer_bucket'],
    output:
        "wes_tx.analysis.txt"
    shell:
        """gsutil -m cp -r analysis/ {params.transfer_bucket} && 
        touch {output}"""

rule transfer_benchmarks:
    input:
        "benchmarks/all_wes_targets.txt" #stub file b/c snkmk doesn't allow dir
    params:
        transfer_bucket=config['transfer_bucket'],
    output:
        "wes_tx.benchmarks.txt"
    shell:
        """gsutil -m cp -r benchmarks/ {params.transfer_bucket} && 
        touch {output}"""

rule transfer_src:
    input:
        "cidc_wes/wes.snakefile" #stub file b/c snkmk doesn't allow dirs
    params:
        transfer_bucket=config['transfer_bucket'],
    output:
        "wes_tx.src.txt"
    shell:
        """gsutil -m cp -r cidc_wes/ {params.transfer_bucket} && 
        touch {output}"""

rule transfer_config_meta:
    input:
        conf="config.yaml",
        meta="metasheet.csv",
    params:
        transfer_bucket=config['transfer_bucket']
    output:
        "wes_tx.config_meta.txt"
    shell:
        """gsutil -m cp -r {input.conf} {params.transfer_bucket} && 
        gsutil -m cp -r {input.meta} {params.transfer_bucket} && 
        touch {output}"""

rule transfer_nohups:
    input:
        "nohup.out"
    params:
        transfer_bucket=config['transfer_bucket'],
        files=lambda wildcards: glob.glob('nohup*.out')
    output:
        "wes_tx.nohup.txt"
    run:
        for f in params.files:
            shell("gsutil -m cp -r {f} {params.transfer_bucket}")
        shell("touch {output}")

