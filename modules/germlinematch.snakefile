def germlinevcf_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append("analysis/germline/%s/%s_SNP92.recode.vcf" % (sample_name, sample_name))
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp

def germ_getNormal(wildcards):
    return germlinevcf_runsHelper(wildcards, 0)

def germ_getTumor(wildcards):
    return germlinevcf_runsHelper(wildcards, 1)

def checkmatch_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        #ls.append("aggregated/%s/a.txt" % (run)),
        ls.append("analysis/germline/%s/%s_matchinformation.txt" % (run,run))
        #ls.append("analysis/germline/%s/matchfile.txt" % run)
    return ls

rule germlinematch_all:  #aggregate of all rule
    input:
        checkmatch_targets

rule vcftoolsdiff:
    input:
        tumor = germ_getTumor,
        normal = germ_getNormal,
    output:
        out="analysis/germline/{run}/{run}.diff.sites_in_files",
	log="analysis/germline/{run}/{run}.log",
    params: run_name=lambda wildcards:wildcards.run, 
    shell:
        "vcftools --vcf  {input.tumor} --diff  {input.normal} --diff-site --out analysis/germline/{params.run_name}/{params.run_name}"

rule testmatch:
    input: 
    	"analysis/germline/{run}/{run}.log"
    output:
        "analysis/germline/{run}/{run}_matchinformation.txt"
    shell:
        "cidc_wes/modules/scripts/Match.py -i {input} > {output}"


