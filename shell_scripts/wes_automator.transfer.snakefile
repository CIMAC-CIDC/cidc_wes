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
    #input:
    #    "analysis/align/mapping.csv" #stub file b/c snkmk doesn't allow dirs
    params:
        transfer_path=config['transfer_path'],
    output:
        "wes_tx.analysis.txt"
    shell:
        """gsutil -m cp -r analysis/ {params.transfer_path} && 
        touch {output}"""

rule transfer_benchmarks:
    #input:
    #    "benchmarks/all_wes_targets.txt" #stub file b/c snkmk doesn't do dir
    params:
        transfer_path=config['transfer_path'],
    output:
        "wes_tx.benchmarks.txt"
    shell:
        """gsutil -m cp -r benchmarks/ {params.transfer_path} && 
        touch {output}"""

rule transfer_src:
    #input:
    #    "cidc_wes/wes.snakefile" #stub file b/c snkmk doesn't allow dirs
    params:
        transfer_path=config['transfer_path'],
    output:
        "wes_tx.src.txt"
    shell:
        """gsutil -m cp -r cidc_wes/ {params.transfer_path} && 
        touch {output}"""

rule transfer_config_meta:
    input:
        conf="config.yaml",
        meta="metasheet.csv",
    params:
        transfer_path=config['transfer_path']
    output:
        "wes_tx.config_meta.txt"
    shell:
        """gsutil -m cp -r {input.conf} {params.transfer_path} && 
        gsutil -m cp -r {input.meta} {params.transfer_path} && 
        touch {output}"""

rule transfer_nohups:
    input:
        "nohup.out"
    params:
        transfer_path=config['transfer_path'],
        files=lambda wildcards: glob.glob('nohup*.out*')
    output:
        "wes_tx.nohup.txt"
    run:
        for f in params.files:
            shell("gsutil -m cp -r {f} {params.transfer_path}")
        shell("touch {output}")

