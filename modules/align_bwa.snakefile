#MODULE: Align fastq files to genome - BWA specific calls
#PARAMETERS:
_logfile="analysis/logs/align.log"
_bwa_threads=16

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

rule bwa_mem:
    input:
        getFastq
    output:
        temp("analysis/align/{sample}/{sample}.bam")
    params:
        index=config['bwa_index'],
        read_group= lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample),
        input_bases="10000000",
    threads: _bwa_threads
    message: "ALIGN: Running BWA mem for alignment"
    log: _logfile
    shell:
    	"""bwa mem -M -R \"{params.read_group}\" -t {threads} -K {params.input_bases} {params.index} {input} | samtools view -Sb - > {output}"""
	#MOVING to jingxin's call
        #"bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"



