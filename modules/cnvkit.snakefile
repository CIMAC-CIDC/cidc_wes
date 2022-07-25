# module: call CNVs with cnvkit
_cnvkit_threads=8

center_cnvkit_ref={"mda": "./ref_files/hg38/cnvkit/cidc_hg38_5k_mdaS1400i_S1609_cnvkit_reference.cnn",
                   "broad": "./ref_files/hg38/cnvkit/cidc_hg38_5k_broad10021_10026ArmB_cnvkit_reference.cnn",
                   "flat": "./ref_files/hg38/cnvkit/cidc_hg38_5k_cnvkitFlat_reference.cnn",
}

def cnvkit_targets(wildcards):
    """ Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        r = config['runs'][run]
        tmr = r[1]
        ls.append("analysis/cnvkit/%s/%s_recalibrated.cns" % (run,tmr))
        ls.append("analysis/cnvkit/%s/%s_recalibrated.call.cns" % (run,tmr))
        ls.append("analysis/cnvkit/%s/%s_recalibrated.call.enhanced.cns" % (run,tmr))
        ls.append("analysis/cnvkit/%s/%s_cnvkit_gainLoss.bed" % (run,run))

        #renamed files for ingestion --see rule cnvkit_rename
        ls.append("analysis/cnvkit/%s/%s.call.cns" % (run,run))
        ls.append("analysis/cnvkit/%s/%s.call.enhanced.cns" % (run,run))
        ls.append("analysis/cnvkit/%s/%s.scatter.png" % (run,run))
        
    return ls

rule cnvkit_all:
    input:
        cnvkit_targets
    benchmark: "benchmarks/cnvkit/cnvkit_all.txt"

def cnvkitInputFn(wildcards):
    #run = config['runs'][wildcards.run]
    #tmr = run[1]
    tmr = wildcards.tmr
    #return tmr's recalibrated bam file
    ls = ["analysis/align/%s/%s_recalibrated.bam" % (tmr, tmr)]
    return ls

def getCNVkit_ref(wildcards):
    #if no cimac_center is defined, default to flat
    cimac = config.get('cimac_center', 'flat')
    if cimac in center_cnvkit_ref:
        return center_cnvkit_ref[cimac]
    else:
        #default to flat ref
        return center_cnvkit_ref['flat']
    
rule cnvkit:
    """call cnvs using cnvkit pipeline"""
    input:
        cnvkitInputFn
    output:
        cns = "analysis/cnvkit/{run}/{tmr}_recalibrated.cns", #cnvkit's first pass calls
        calls_cns = "analysis/cnvkit/{run}/{tmr}_recalibrated.call.cns", #cnvkit's more refined calls
        scatter = "analysis/cnvkit/{run}/{tmr}_recalibrated-scatter.png",
    threads: _cnvkit_threads
    group: "cnvkit"
    params:
        cnvkit_ref = lambda wildcards: getCNVkit_ref(wildcards),
        output_dir=lambda wildcards: "analysis/cnvkit/%s" % wildcards.run,
    log: "analysis/logs/cnvkit/{run}/{tmr}.cnvkit.log"
    benchmark:
        "benchmarks/cnvkit/{run}/{tmr}.cnvkit.txt"
    shell:
        """cnvkit.py batch {input} -r {params.cnvkit_ref} -p {threads} --scatter --diagram -d {params.output_dir}"""


def purity_checker(run):
    tumor = config['runs'][run][1]
    normal = config['runs'][run][0]
    output = "analysis/cnvkit/%s/%s_recalibrated.call.enhanced.cns" % (run,tumor)
    vcf="analysis/somatic/%s/%s_tnscope.output.vcf.gz" % (run,run)
    cns="analysis/cnvkit/%s/%s_recalibrated.call.cns" % (run,tumor)
    if 'purity' not in config['skipped_modules']:
        file="analysis/purity/%s/%s.optimalpurityvalue.txt" % (run,run)
        if os.path.exists(file):
            df = pd.read_csv(file, na_filter=False, delimiter="\t")
            if df["purity"][0] != "NA":
                with open(file) as f:
                    tmp = f.readline()
                    purity = "-m clonal --purity %s" % f.readline().split("\t")[2]
            else:
                purity = ""
    return purity #either "-m clonal --purity VAL" or ""

rule cnvkit_enhance:
    """Add somatic snp and purity information to cnvkit's refined call"""
    input:
        cns="analysis/cnvkit/{run}/{tmr}_recalibrated.call.cns",
        vcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz",
        purity="analysis/purity/{run}/{run}.optimalpurityvalue.txt",
    output:
        "analysis/cnvkit/{run}/{tmr}_recalibrated.call.enhanced.cns"
    group: "cnvkit"
    params:
        tmr_name = lambda wildcards:"-i %s" % config['runs'][wildcards.run][1],
	nrm_name = lambda wildcards:"-n %s" % config['runs'][wildcards.run][0] if not config.get('tumor_only') else "",
        run = lambda wildcards: wildcards.run,
        purity=lambda wildcards: purity_checker(wildcards.run)
    log: "analysis/logs/cnvkit/{run}/{tmr}.cnvkit_enhance.log"
    benchmark:
        "benchmarks/cnvkit/{run}/{tmr}.cnvkit_enhance.txt"
    shell:
        "cnvkit.py call {input.cns} -y -v {input.vcf} {params.tmr_name} {params.nrm_name} {params.purity} -o {output}"



# #if 'purity' not in config['skipped_modules']: #run cnvkit call w/ purity
# if purity_checker(wildcards): #run cnvkit call w/ purity
#     rule cnvkit_enhance:
#         """Add somatic snp and purity information to cnvkit's refined call"""
#         input:
#             cns="analysis/cnvkit/{run}/{tmr}_recalibrated.call.cns",
#             vcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz",
#             purity="analysis/purity/{run}/{run}.optimalpurityvalue.txt",
#         output:
#             "analysis/cnvkit/{run}/{tmr}_recalibrated.call.enhanced.cns"
#         group: "cnvkit"
#         params:
#             tmr_name = lambda wildcards:"-i %s" % config['runs'][wildcards.run][1],
#             #check for tumor-only?
#             nrm_name = lambda wildcards:"-n %s" % config['runs'][wildcards.run][0] if not config.get('tumor_only') else "",
#         log: "analysis/logs/cnvkit/{run}/{tmr}.cnvkit_enhance.log"
#         benchmark:
#             "benchmarks/cnvkit/{run}/{tmr}.cnvkit_enhance.txt"
#         shell:
#             #first cmd grabs the purity value, 3rd col of 2nd line
#             """PURITY=$(sed -n 2p {input.purity} | cut -f 3) && cnvkit.py call {input.cns} -y -v {input.vcf} {params.tmr_name} {params.nrm_name} -m clonal --purity $PURITY -o {output}"""

# else: #run cnvkit call w/o purity

#     rule cnvkit_enhance_noPurity:
#         """Add somatic snp and purity information to cnvkit's refined call"""
#         input:
#             cns="analysis/cnvkit/{run}/{tmr}_recalibrated.call.cns",
#             vcf="analysis/somatic/{run}/{run}_tnscope.output.vcf.gz",
#         output:
#             "analysis/cnvkit/{run}/{tmr}_recalibrated.call.enhanced.cns"
#         group: "cnvkit"
#         params:
#             tmr_name = lambda wildcards:"-i %s" % config['runs'][wildcards.run][1],
#             #check for tumor-only?
#             nrm_name = lambda wildcards:"-n %s" % config['runs'][wildcards.run][0] if not config.get('tumor_only') else "",
#         log: "analysis/logs/cnvkit/{run}/{tmr}.cnvkit_enhance_noPurity.log"
#         benchmark:
#             "benchmarks/cnvkit/{run}/{tmr}.cnvkit_enhance_noPurity.txt"
#         shell:
#             """cnvkit.py call {input.cns} -y -v {input.vcf} {params.tmr_name} {params.nrm_name} -o {output}"""

def cnvkit_callGainLossInput(wildcards):
    run = config['runs'][wildcards.run]
    tmr = run[1]
    ls = ["analysis/cnvkit/%s/%s_recalibrated.call.enhanced.cns" % (wildcards.run, tmr)]
    return ls

rule cnvkit_callGainLoss:
    """use hard-cutoffs to call regions of GAIN/LOSS"""
    input:
        cnvkit_callGainLossInput
    output:
        #NOTE: changing from {run}-{tmr} to {run}-{run} to be more consistent
        "analysis/cnvkit/{run}/{run}_cnvkit_gainLoss.bed"
    group: "cnvkit"
    log: "analysis/logs/cnvkit/{run}/{run}.cnvkit_callGainLoss.log"
    benchmark:
        "benchmarks/cnvkit/{run}/{run}.cnvkit_callGainLoss.txt"
    shell:
        "./cidc_wes/modules/scripts/copynumber_callGainLoss.py -f {input} -o {output}"

def cnvkit_renameInput(wildcards):    
    run = wildcards.run
    tmr = config['runs'][run][1]
    
    cns="analysis/cnvkit/%s/%s_recalibrated.call.cns" % (run, tmr)
    enh_cns="analysis/cnvkit/%s/%s_recalibrated.call.enhanced.cns" % (run, tmr)
    scatter="analysis/cnvkit/%s/%s_recalibrated-scatter.png" % (run, tmr)
    tmp = {'cns': cns, 'enhanced_cns': enh_cns, 'scatter': scatter}
    return tmp

#RENAME the output to canonical {run}/{run} format so that ingestion is easier
rule cnvkit_rename:
    input:
        unpack(cnvkit_renameInput)
    output:
        cns="analysis/cnvkit/{run}/{run}.call.cns",
        enhanced_cns="analysis/cnvkit/{run}/{run}.call.enhanced.cns",
        scatter="analysis/cnvkit/{run}/{run}.scatter.png",
    shell:
        """cp {input.cns} {output.cns} && cp {input.enhanced_cns} {output.enhanced_cns} && cp {input.scatter} {output.scatter}"""
