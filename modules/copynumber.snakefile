#MODULE: CNV  calls by Sentieon

_cnvcall_threads=32

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

def cnv_getNormal_sample(wildcards):
    return cnvcall_runsHelper(wildcards, 0)

def cnv_getTumor_sample(wildcards):
    return cnvcall_runsHelper(wildcards, 1)

def copynumber_getCohortFiles(wildcards):
    """Returns the analysis/copynumber/{run}/{run} for a given cohort"""
    c = wildcards.cohort
    ls = []
    for run in config['cohorts'][c]:
        ls.append("analysis/copynumber/%s/%s_cnvcalls.txt" % (run,run))
    return ls

def copynumber_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #ls.append("analysis/copynumber/%s/%s_pon.hdf5" % (run,run))
        ls.append("analysis/copynumber/%s/%s_cnvcalls.txt"%(run,run))
        ls.append("analysis/copynumber/%s/%s_cnvcalls.circos.txt"%(run,run))

    #Check for cohorts
    if config['cohorts']:
        for cohort in config['cohorts']:
            ls.append("analysis/cohorts/%s/%s_merged.cnvcalls.txt" % (cohort,cohort))
            #Gistic output here
            #ls.append("analysis/cohorts/%s/%s_merged.cnvcalls.txt" % (cohort,cohort))
    return ls

rule copynumber_all:
    input:
        copynumber_targets

#NOT NEEDED?--leave in for now just in case we want users to be able to do this
rule copynumber_create_pon_sentieon:
    input:
        normal_recalibratedbam="analysis/align/{sample}/{sample}_recalibrated.bam" ##get only Normal samples
    output:
        ponfile="analysis/copynumber/{run}/{run}_pon.hdf5"
    params:
        index=config['genome_fasta'],
        index1=config['sentieon_path'],
        #target bed provided by user- not enabled for now
        #targetbed= config['pon_target_bed'],
    threads: _cnvcall_threads
    benchmark:
        "benchmarks/CNV/{run}/{run}.create_pon_sentieon.txt"
    shell:
        """{params.index1}/sentieon driver  -t {threads} -r {params.index} -i {input.normal_recalibratedbam} --algo CNV  --target {input.targetbed} --target_padding 0 --create_pon {output.ponfile}"""


rule copynumber_CNVcall:
    input:
        #ONLY perform this analysis for Tumor samples-
        tumor_recalibratedbam = cnv_getTumor_sample
    output:
        cnvcalls="analysis/copynumber/{run}/{sample}_cnvcalls.txt",
        tWeights="analysis/copynumber/{run}/{sample}_cnvcalls.txt.targetWeights",
        tsv="analysis/copynumber/{run}/{sample}_cnvcalls.txt.tn.tsv",
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        ponfile=config['pons'],
        target=config['pons_target'],
    group: "copynumber"
    threads:_cnvcall_threads
    benchmark:
        "benchmarks/copynumber/{run}/{sample}.CNVcall.txt"
    shell:
        #NOTE: target_padding is being set to 0--following the broad's method
        #"""{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.tumor_recalibratedbam} --algo CNV --target {params.target} --target_padding 0  --pon {params.ponfile} {output.cnvcalls}"""
        """{params.sentieon_path}/sentieon driver -t {threads} -r {params.index} -i {input.tumor_recalibratedbam} --algo CNV  --pon {params.ponfile} {output.cnvcalls}"""

rule copynumber_processCNVcircos:
    """Process the cnv output file to be an adaquate input into the circos 
    plot"""
    input:
        cnvcalls="analysis/copynumber/{run}/{sample}_cnvcalls.txt.tn.tsv",
    output:
        cnvcalls="analysis/copynumber/{run}/{sample}_cnvcalls.circos.txt",
    group: "copynumber"
    benchmark:
        "benchmarks/copynumber/{run}/{sample}.processCNVcircos.txt"
    shell:
        "cidc_wes/modules/scripts/copynumber_processCircos.py -f {input} > {output}"

rule copynumber_mergeCohorts:
    """Merge the cnvcalls.txt files from the runs in a cohort"""
    input:
        copynumber_getCohortFiles
    params:
        files = lambda wildcards,input: " -f ".join(input)
    group: "copytnumber"
    output:
        "analysis/cohorts/{cohort}/{cohort}_merged.cnvcalls.txt"
    benchmark:
        "benchmarks/cohorts/{cohort}/{cohort}_copynumber_mergeCohorts"
    shell:
        "cidc_wes/modules/scripts/copynumber_mergeCohorts.py -f {params.files} > {output}"

rule copynumber_gistic2:
    """Call gistic2 to convert from segment-based cnv calls to gene-based"""
    input:
        "analysis/cohorts/{cohort}/{cohort}_merged.cnvcalls.txt"
    params:
        #ref: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/
        ref=config['gistic2'],
        ta= 0.1, #amplification threshold
        td=0.1, #deletion threshold
        armpeel=1, #perform arm level peel off
        brlen=0.7, #broad.length cutoff
        cap=1.5, #Not in gistic2 doc
        conf=0.99, #confidence level
        genegistic=1, #turn on gene level for deletions
        gcm="extreme", #not in gistic2 doc
        js=4, #join segment
        maxseg=2000, #max.sample.segment
        qvt=0.25, #qv threshold
        rx=0, #remove x
        savegene=1, #not in gistic2 doc
        #HARD-code gistic2 paths--need a better solution
        ld_lib_path="export LD_LIBRARY_PATH=/usr/local/bin/gistic2/MCR_Installer/v83/runtime/glnxa64:/usr/local/bin/gistic2/MCR_Installer/v83/bin/glnxa64:/usr/local/bin/gistic2/MCR_Installer/v83/sys/os/glnxa64:",
        xapp_dir="export XAPPLRESDIR=/usr/local/bin/gistic2/MCR_Installer/v83/X11/app-defaults",

        outdir =lambda wildcards,input,output: "/".join(str(output).split("/")[:-1])
    group: "copytnumber"
    output:
        "analysis/cohorts/{cohort}/gistic2/all_data_by_genes.txt"
    benchmark:
        "benchmarks/cohorts/{cohort}/{cohort}_copynumber_gistic2"
    shell:
        "{params.ld_lib_path};{params.xapp_dir}; gistic2-bin -refgene {params.ref} -seg {input} -b {params.outdir} -ta {params.ta} -armpeel {params.armpeel} -brlen {params.brlen} -cap {params.cap} -conf {params.conf} -td {params.td} -genegistic {params.genegistic} -gcm {params.gcm} -js {params.js} -maxseg {params.maxseg} -qvt {params.qvt} -rx {params.rx} -savegene {params.savegene}"
