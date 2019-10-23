#module: neoantigen prediction module using pvacseq pipeline (from pvactools)

_neoantigen_threads=64 #should be set to as max cores; w/ 64 runtime~=1hr

def neoantigen_runsHelper(wildcards, iindex):
    """Given a snakemake wildcards, an iindex - 0 for Normal, 1 for Tumor,
    returns the sample name of Normal (if iindex=0) else sample name of Tmr"""
    tmp = []
    r = config['runs'][wildcards.run]
    #print(r)

    #check that we have a valid pair
    if len(r) >=2:
        sample_name = r[iindex]
        tmp.append(sample_name) 
    else:
        #NOTE: I can't figure out a proper kill command so I'll do this
        tmp=["ERROR! BAD pairing for run--requires at least two samples: %s" % (wildcards.run)]
    #print(tmp)
    return tmp


def neoantigen_getNormal(wildcards):
    return neoantigen_runsHelper(wildcards, 0)

def neoantigen_getTumor(wildcards):
    return neoantigen_runsHelper(wildcards, 1)

def neoantigen_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config['runs']:
        tumor = config['runs'][run][1] 
        ls.append("analysis/neoantigen/%s/MHC_Class_I/%s.filtered.condensed.ranked.tsv" % (run,tumor))
        ls.append("analysis/neoantigen/%s/MHC_Class_I/%s.filtered.condensed.ranked.addSample.tsv" % (run,tumor))
    return ls

def getPvacseqOut(wildcards):
    """returns tuple (run, tumor sample name)"""
    run = wildcards.run
    tumor = config['runs'][run][1]
    return(run,tumor)

def getTumorHLA(wildcards):
    """get the optitype results file for the tumor sample"""
    run = wildcards.run
    tumor = config['runs'][run][1]
    ls = ["analysis/optitype/%s/%s_result.tsv" % (tumor, tumor)]
    if 'neoantigen_run_classII' in config and config['neoantigen_run_classII']:
        ls.append("analysis/xhla/%s/report-%s-hla.json" % (tumor,tumor))
    #print(ls)
    return ls

def parseHLA(hla_files):
    """Given an optitypes '_results.tsv' file; parses the HLA A, B, C
    and returns these as a comma-separated string (for pvacseq) input

    NOTE: cureently the optitype results.tsv looks somthing like this:
    	A1	A2	B1	B2	C1	C2	Reads	Objective
    0					C*06:04	C*06:04	4.0	3.99
    **So were' going to parse cols 1-6 and return that"""

    #CATCH when the HLA does not exist yet
    #print(optitype_out_file)
    optitype_out_file = hla_files[0]
    if not os.path.exists(optitype_out_file):
        #print("WES WARNING: %s is not found!" % optitype_out_file)
        return ""

    f = open(optitype_out_file)
    hdr = f.readline().strip().split("\t") #ignore for now
    classI = f.readline().strip().split("\t")[1:7] #want first 6 cols
    #FOR classI alleles, prepend a HLA to each of them
    classI = ["HLA-%s" % a for a in classI if a]
    #print(classI)
    f.close()
    
    #check for xhla file
    classII = []
    if 'neoantigen_run_classII' in config and config['neoantigen_run_classII'] and len(hla_files) > 1:
        xhla_out_file = hla_files[1]
        
        #PARSE xhla json file...
        f = open(xhla_out_file)
        xhla_out = json.load(f)
        f.close()

        #build classII alleleles
        #ONLY add class II alleles--i.e. ones that start with "D"
        classII = [a for a in xhla_out['hla']['alleles'] if a.startswith("D")]
        #print(classII)
        
    if classII:
        classI.extend(classII)
    #NOTE: NOW classI has all hla alleles (including classII if opted for)
    hla = ",".join(["%s" % a for a in classI if a])
    #print(hla)
    return hla
    
    

rule neoantigen_all:
    input:
        neoantigen_targets

rule neoantigen_vep_annotate:
    input:
        "analysis/somatic/{run}/{run}_tnsnv.filter.vcf"
    output:
        "analysis/somatic/{run}/{run}_tnsnv.filter.neoantigen.vep.vcf"
    params:
        index=config['genome_fasta'],
        vep_data=config['vep_data'],
        vep_plugins=config['vep_plugins'],

        #normal = lambda wildcards: config['runs'][wildcards.run][0],
        #tumor = lambda wildcards: config['runs'][wildcards.run][1],
    group: "neoantigen"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/neoantigen/{run}/{run}.neoantigen_vep_annotate.txt"
    shell:
        """vep --input_file {input} --output_file {output} --format vcf --vcf --symbol --terms SO --tsl --hgvs --fasta {params.index} --offline --cache --dir_cache {params.vep_data} --plugin Downstream --plugin Wildtype --dir_plugins {params.vep_plugins} --pick"""


rule neoantigen_pvacseq:
    """NOTE: neoantigen's pvacseq is not available on CONDA
    MUST either be installed in base system/docker container"""
    input:
        vcf="analysis/somatic/{run}/{run}_tnsnv.filter.neoantigen.vep.vcf",
        hla=getTumorHLA,
    output:
        main="analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.condensed.ranked.tsv",
        #OTHERS:
        filtered="analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.tsv",
        all_epitopes="analysis/neoantigen/{run}/MHC_Class_I/{tumor}.all_epitopes.tsv",
        #The output below causes a rule ambiguity with neoantigen_addSample's
        #output file
        #tsv="analysis/neoantigen/{run}/MHC_Class_I/{tumor}.tsv",
    
        #NOTE: can't get this last one in b/c snakemake complains about
        #a non-unique name
        #input_log="analysis/neoantigen/{run}/MHC_Class_I/log/inputs.yml",
    
        #NOTE: the wildcard tumor is used b/c we can't actually pull in the
        #tumor name 

        #NOTE: typically when HLA class I and HLA class II are both called
        #the results are in combined/{inputname}.condensed.tsv
        #Since we're only generating class I for now, no combined is generated
    params:
        normal = lambda wildcards: config['runs'][wildcards.run][0],
        tumor = lambda wildcards: config['runs'][wildcards.run][1],
        iedb = config['neoantigen_iedb'],
        HLA = lambda wildcards,input: parseHLA(input.hla),
        callers=config['neoantigen_callers'] if config['neoantigen_callers'] else "MHCflurry NetMHCcons MHCnuggetsII",
        epitope_lengths=config['neoantigen_epitope_lengths'] if config['neoantigen_epitope_lengths'] else "8,9,10,11",
        output_dir = lambda wildcards: "%sanalysis/neoantigen/%s/" % (config['remote_path'], wildcards.run),
    threads: _neoantigen_threads
    group: "neoantigen"
    log: "analysis/logs/neoantigen/{run}/{tumor}.neoantigen_pvacseq.log"
    benchmark:
        "benchmarks/neoantigen/{run}/{tumor}.neoantigen_pvacseq.txt"
    shell:
        """pvacseq run {input.vcf} {params.tumor} {params.HLA} {params.callers} {params.output_dir} -e {params.epitope_lengths} -t {threads} --normal-sample-name {params.normal} --iedb-install-directory {params.iedb} 2> {log}"""

  
rule neoantigen_add_sample:
    """Takes the {tumor}.filtered.condensed.ranked.tsv file and adds a first
    col, Sample which is set to the RUN name (not the sample name!)
    to be used for pvacseq_plot.R"""
    input:
        "analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.condensed.ranked.tsv"
    output:
        "analysis/neoantigen/{run}/MHC_Class_I/{tumor}.filtered.condensed.ranked.addSample.tsv"
    params:
        run_name = lambda wildcards: wildcards.run
    group: "neoantigen"
    shell:
        "cidc_wes/modules/scripts/neoantigen_addSample.py -f {input} -n {params.run_name} -o {output}"

