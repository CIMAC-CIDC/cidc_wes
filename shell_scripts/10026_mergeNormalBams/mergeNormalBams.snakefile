##10026 Trial Merge Patient and Donor Germline bam files
_threads = 4
configfile: "config.yaml"


#generates a dictionary of sub-files -> merged names to help us when we
#need to reheader the sub-files to have the merged file names
def mergedSampleNameMap(config):
    tmp = {}
    for sample in config["samples"]:
        for subBam in config["samples"][sample]:
            tmp[subBam] = sample
    return tmp
_nameMap = mergedSampleNameMap(config)
    
def targets(wildcards):
    ls = []
    for sample in config["samples"]:
        ls.append("analysis/merged/%s.bam" % sample)
    return ls

rule all:
    input:
      targets


def getSubBam(wildcards):
    subBam = wildcards.subBam
    sample = _nameMap[subBam]
    bam = config['samples'][sample][subBam]
    #print(bam)
    return bam

rule newHeader_subBam:
    """Replaces the subBam name with the merged file name in the bam header
    outputs new header
    """
    input:
        getSubBam
    output:
        "analysis/headers/{subBam}.header.txt"
    params:
        #replace subBam header string with new merged file name
        sed_cmd=lambda wildcards: "\'s/SM:%s/SM:%s/g\'" % (wildcards.subBam, _nameMap[wildcards.subBam]),
    shell:
        #Get the header and replace the sub bam sample name in @RG tags
        """samtools view -H {input} | sed {params.sed_cmd} > {output}"""

def reheader_subBam_input(wildcards):
    subBam = wildcards.subBam
    sample = _nameMap[subBam]
    bam = config['samples'][sample][subBam]
    header = "analysis/headers/%s.header.txt" % subBam
    tmp = {'bam': bam, 'header': header}
    return tmp

rule reheader_subBam:
    input:
        unpack(reheader_subBam_input)
    output:
        temp("analysis/reheader/{subBam}.bam")
    shell:
        "samtools reheader {input.header} {input.bam} > {output}"

def merge_subBams_input(wildcards):
    sample = wildcards.sample
    ls = []
    for subBam in config['samples'][sample]:
        ls.append("analysis/reheader/%s.bam" % subBam)
    return ls

rule merge_subBams:
    input:
        merge_subBams_input
    output:
        "analysis/merged/{sample}.bam"
    threads: _threads
    shell:
        """sambamba merge -t {threads} {output} {input}"""
