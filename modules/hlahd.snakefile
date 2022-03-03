# module: HLA-HD HLA caller
_hlahd_threads=8

def hlahd_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for sample in config['samples']:
        ls.append("analysis/hlahd/%s/result/%s_final.result.txt" % (sample, sample))
    return ls

rule hlahd_all:
    input:
        hlahd_targets
    benchmark: "benchmarks/hlahd/hlahd_all.txt"

rule hlahd:
    """calculate hlatyping by hla-hd"""
    input:
        chr6fastqfile1 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end1.fastq",
        chr6fastqfile2 = "analysis/optitype/{sample}/{sample}.sorted.chr6.end2.fastq"
    output:
        "analysis/hlahd/{sample}/result/{sample}_final.result.txt"
    threads: _hlahd_threads
    group: "hlahd"
    params:
        name=lambda wildcards: wildcards.sample,
        output_dir=lambda wildcards: "%sanalysis/hlahd/" % config['remote_path'],
        min_lengh = 50, #hla-hd param -m which sets the min read length- fixed b/c we don't expect any libraries with shorter than 50bp
        cut_perc = 0.95, #hla-hd param -c, if mistmatch how much to cut back
        freq_data = config['hlahd_freq_data'],
        split_file = config['hlahd_split_file'],
        dictionary = config['hlahd_dictionary'],
    log: "analysis/logs/hlahd/{sample}/{sample}.hlahd.log"
    benchmark:
        "benchmarks/hlahd/{sample}/{sample}.hlahd.txt"
    shell:
        """hlahd.sh -m {params.min_lengh} -c {params.cut_perc} -t {threads} -f {params.freq_data} {input.chr6fastqfile1} {input.chr6fastqfile2} {params.split_file} {params.dictionary} {params.name} {params.output_dir} 2> {log}"""
