#MODULE: CNV  calls by Sentieon

_cnvcall_threads=8

def cnvcall_runsHelper(wildcards,pairs):
    """Given a snakemake wildcards,  pairs - 0 for Normals, 1 for Tumors,
    returns the  Normal (if pairs=0) else Tumor"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[pairs]
        tmp.append("analysis/align/%s/%s_recalibrated.bam" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp

def getNormal_sample(wildcards):
    return cnvcall_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return cnvcall_runsHelper(wildcards, 1)

def cnvcall_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #ls.append("analysis/cnvcalls/%s/%s_pon.hdf5" % (run,run))
        ls.append("analysis/cnvcalls/%s/%s_cnvcalls.txt"%(run,run))
    return ls

rule cnvcalls_all:
    input:
        cnvcall_targets

#NOT NEEDED?--leave in for now just in case we want users to be able to do this
rule create_pon_sentieon:
    input:
        normal_recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam" ##get only Normal samples
    output:
        ponfile="analysis/cnvcalls/{run}/{run}_pon.hdf5"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        #target bed provided by user- not enabled for now
        #targetbed= config['pon_target_bed'],
    threads: 32 #_cnvcall_threads
    benchmark:
        "benchmarks/CNV/{sample}/{sample}.create_pon_sentieon.txt"
    shell:
        """{params.index1}/sentieon driver  -t {threads} -r {params.index} -i {input.normal_recalibratedbam} --algo CNV  --target {input.targetbed} --target_padding 0 --create_pon {output.ponfile}"""


rule CNVcall_sentieon:
    input:
        #ONLY perform this analysis for Tumor samples-
        tumor_recalibratedbam = getTumor_sample
    output:
        cnvcalls="analysis/cnvcalls/{run}/{sample}_cnvcalls.txt"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        ponfile=config['pons'],
        target=config['pons_target'],
    threads:_cnvcall_threads
    benchmark:
        "benchmarks/CNV/{sample}/{sample}.cnvcall.txt"
    shell:
        #NOTE: target_padding is being set to 0--following the broad's method
        #"""{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.tumor_recalibratedbam} --algo CNV --target {params.target} --target_padding 0  --pon {params.ponfile} {output.cnvcalls}"""
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.tumor_recalibratedbam} --algo CNV  --pon {params.ponfile} {output.cnvcalls}"""
    
