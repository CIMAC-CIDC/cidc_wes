#MODULE: CNV  calls by Sentieon


_cnvcall_threads=64


def cnvcall_runsHelper(wildcards,pairs):
    """Given a snakemake wildcards,  pairs - 0 for Normals, 1 for Tumors,
    returns the  Normal (if pairs=0) else Tumor"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[pairs]
        tmp.append(sample_name)
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def getNormal_sample(wildcards):
    return somatic_runsHelper(wildcards, 0)

def getTumor_sample(wildcards):
    return somatic_runsHelper(wildcards, 1)

def cnvcall_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #ONLY perform this analysis for Tumor samples-
        tumor = getTumor_sample(wildcards)
        #ls.append("analysis/cnvcalls/%s/%s_pon.hdf5" % (run,run))
        ls.append("analysis/cnvcalls/%s/%s_cnvcalls"%(tumor,tumor))
    return ls

rule cnvcalls_all:
    input:
        cnvcall_targets

#NOT NEEDED?
# rule create_pon_sentieon:
#     input:
#         #target bed provided by broad
#         targetbed="/cluster/jxfu/proj/CIDC/Sentieon/data/Reference/broad_revised.bed",
#         normal_recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam" ##get only Normal samples
#     output:
#         ponfile="analysis/cnvcalls/{run}/{run}_pon.hdf5"
#     params:
#         index=config['genome_fasta'],
#         index1=config['sentieon_path']
#     threads:_cnvcall_threads
#     shell:
#         """{params.index1}/sentieon driver  -t {threads} -r {params.index} -i {input.normal_recalibratedbam} --algo CNV  --target {input.targetbed} --target_padding 0 --create_pon {output.ponfile}"""


rule CNVcall_sentieon:
    input:
        tumor_recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam",        
    output:
        cnvcalls="analysis/cnvcalls/{sample}/{sample}_cnvcalls"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        ponfile=config['PoNs'],
        target=config['PoNs_target'],
    threads:_cnvcall_threads
    shell:
        #NOTE: target_padding is being set to 0--is this right?
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.tumor_recalibratedbam} --algo CNV --target {params.target} --target_padding 0  --pon {params.ponfile} {output.cnvcalls}"""
    
