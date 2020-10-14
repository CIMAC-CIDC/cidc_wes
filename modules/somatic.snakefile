#module: Somatic Variant calls by Sentieon
#import os
#from string import Template

_lego_plot_data_path="cidc_wes/cidc-vs/cidcvs/data/REF/"
#MAP of TCGA Cancer type names and their data files
_lego_plot_map={'ACC':'TCGA-ACC.mtrc',
                'BLCA':'TCGA-BLCA.mtrx',
                'BRCA':'TCGA-BRCA.mtrx',
                'CESC':'TCGA-CESC.mtrx',
                'CHOL':'TCGA-CHOL.mtrx',
                'COAD':'TCGA-COAD.mtrx',
                'DLBC':'TCGA-DLBC.mtrx',
                'ESCA':'TCGA-ESCA.mtrx',
                'GBM':'TCGA-GBM.mtrx',
                'HNSC':'TCGA-HNSC.mtrx',
                'KICH':'TCGA-KICH.mtrx',
                'KIRC':'TCGA-KIRC.mtrx',
                'KIRP':'TCGA-KIRP.mtrx',
                'LGG':'TCGA-LGG.mtrx',
                'LIHC':'TCGA-LIHC.mtrx',
                'LUAD':'TCGA-LUAD.mtrx',
                'LUSC':'TCGA-LUSC.mtrx',
                'MESO':'TCGA-MESO.mtrx',
                'OV':'TCGA-OV.mtrx',
                'PAAD':'TCGA-PAAD.mtrx',
                'PANCAN':'TCGA-PANCAN.mtrx',
                'PCPG':'TCGA-PCPG.mtrx',
                'PRAD':'TCGA-PRAD.mtrx',
                'READ':'TCGA-READ.mtrx',
                'SARC':'TCGA-SARC.mtrx',
                'SKCM':'TCGA-SKCM.mtrx',
                'STAD':'TCGA-STAD.mtrx',
                'TGCT':'TCGA-TGCT.mtrx',
                'THCA':'TCGA-THCA.mtrx',
                'THYM':'TCGA-THYM.mtrx',
                'UCEC':'TCGA-UCEC.mtrx',
                'UCS':'TCGA-UCS.mtrx',
                'UVM':'TCGA-UVM.mtrx'}


_somatic_threads=32
#_vcf2maf_threads=4

def build_tcga_param():
    """Builds the TCGA panel parameter to use for the mutationSignature rule"""
    #PANCAN is a permanent member of the panel
    ls = []
    tcga_panel = config.get('tcga_panel')

    if tcga_panel:
        for cancer in tcga_panel.split(" "):
            if cancer in _lego_plot_map:
                ls.append(os.path.join(_lego_plot_data_path, _lego_plot_map[cancer]))
    ls.append(os.path.join(_lego_plot_data_path, _lego_plot_map['PANCAN']))
    tmp = " -c ".join(ls)
    return tmp #need to add a -c

#NOTE: somatic_runsHelper, getNormal_sample, and getTumor_sample are NOT
#called by any one!
def somatic_runsHelper(wildcards, iindex):
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


def somatic_getNormal(wildcards):
    return somatic_runsHelper(wildcards, 0)

def somatic_getTumor(wildcards):
    return somatic_runsHelper(wildcards, 1)

def somatic_getNormal_recal(wildcards):
    sample = somatic_runsHelper(wildcards, 0)[0]
    return "analysis/align/%s/%s_recalibrated.bam" % (sample,sample)

def somatic_getNTumor_recal(wildcards):
    sample = somatic_runsHelper(wildcards, 1)[0]
    return "analysis/align/%s/%s_recalibrated.bam" % (sample,sample)

def somatic_getNormal_recal_bai(wildcards):
    sample = somatic_runsHelper(wildcards, 0)[0]
    return "analysis/align/%s/%s_recalibrated.bam.bai" % (sample,sample)

def somatic_getNTumor_recal_bai(wildcards):
    sample = somatic_runsHelper(wildcards, 1)[0]
    return "analysis/align/%s/%s_recalibrated.bam.bai" % (sample,sample)

def somatic_helper_targets(wildcards, caller):
    ls = []
    center = config.get('cimac_center', 'broad') #Try to get center, default broad
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s_%s.output.vcf.gz" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.vcf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.vcf.gz" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.exons.vcf.gz" % (run,run, caller))
        #next 3 for filter.maf/pdf
        ls.append("analysis/somatic/%s/%s_%s.filter.vep.vcf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.maf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.pdf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.stats.txt" % (run,run, caller))
        #next 2 for mutation load
        ls.append("analysis/somatic/%s/%s_%s.output.exon.maf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.mutationload.txt" % (run,run, caller))
        #next 2 for circos
        ls.append("analysis/somatic/%s/%s_%s.indel.circos.txt" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.snp.circos.txt" % (run,run, caller))

        ls.append("analysis/somatic/%s/%s_%s_somatic_SNV_summaries.csv" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.exons.%s.vcf.gz" % (run,run, caller, center))
    ls.append("analysis/somatic/somatic_mutation_summaries.%s.csv" % caller)
    ls.append("analysis/somatic/somatic_functional_annot_summaries.%s.csv" % caller)

    #json files
    ls.append("analysis/report/json/somatic/%s_%s.filtered_maf.json" % (run, caller))
    return ls

def somatic_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    caller = config.get('somatic_caller', 'tnscope')
    ls = somatic_helper_targets(wildcards, caller)

    #output file map
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s.somatic.output.yaml" % (run,run))
        
    return ls

def somatic_output_files(wildcards):
    """returns a list of filepaths generated by this module to store 
    in the CIDC for a given sample 
    """
    ls = []
    caller = config['somatic_caller']
    center = config.get('cimac_center', 'broad')
    run = wildcards.run
    ls.append("analysis/somatic/%s/%s_%s.filter.exons.%s.vcf.gz" % (run,run,caller,center))
    ls.append("analysis/somatic/%s/%s_%s.filter.maf" % (run,run,caller))
    ls.append("analysis/somatic/%s/%s_%s.filter.vcf" % (run,run,caller))
    ls.append("analysis/somatic/%s/%s_%s.output.maf" % (run,run,caller))
    ls.append("analysis/somatic/%s/%s_%s.output.vcf" % (run,run,caller))
    
    return ls

def somatic_make_file_map_makeKeys():
    """Makes the keys for the yaml_writer"""
    center = config.get('cimac_center', 'broad')
    tmp = " -k ".join(['%s_filter_exon' % center,'filter_maf','filter_vcf','output_maf','output_vcf']),

rule somatic_make_file_map:
    input:
        somatic_output_files
    output:
        "analysis/somatic/{run}/{run}.somatic.output.yaml"
    params:
        run = lambda wildcards: wildcards.run,
        #keys = " -k ".join(['broad_filter_exons','mda_filter_exons','mocha_filter_exon','filter_maf','filter_vcf','output_maf','output_vcf']),
        keys = lambda wildcards: somatic_make_file_map_makeKeys(),
        files = lambda wildcards, input: " -f ".join(input),
    shell:
        "cidc_wes/modules/scripts/yaml_writer.py -t runs -n {params.run} -k {params.keys} -f {params.files} > {output}"

rule somatic_all:
    input:
        somatic_targets

###############################################################################
# PLEASE LOOK AT the different caller snakefiles,
# somatic_tnhaplotyper2, somatic_tnsnv, somatic_tnscope
# for caller specific rules!
###############################################################################

rule filter_raw_vcf:
    """General rule to filter the three different types of vcf.gz files"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    params:
        index=config['genome_fasta'],
        sentieon_path=config['sentieon_path'],
        tumor=lambda wildcards: somatic_getTumor(wildcards),
        normal= lambda wildcards: somatic_getNormal(wildcards),
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_filter_raw_vcf.txt"
    run: #DISABLES the conda env
        #SWITCH for tnhaplotyper2 filter
        if (wildcards.caller == "tnhaplotyper2"):
            shell("{params.sentieon_path}/sentieon tnhapfilter --tumor_sample {params.tumor} --normal_sample {params.normal} -v {input} {output}")
        else:
            shell("""vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}""")


rule gunzip_vcf:
    """General rule to gunzip the three types of vcf.gz files-
    tnscope_, tnsnv, and tnhaplotyper"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz"
    output:
        #Should we make this a temp?
        "analysis/somatic/{run}/{run}_{caller}.output.vcf"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_gunzip_vcf.txt"
    group: "somatic"
    shell:
        #NOTE: we want to keep the original .gz vcf file
        "gunzip < {input} > {output}"

rule vcfVEP:
    """Rule to annotate vcf files with vep"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.{type}.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.vep.vcf"
    params:
        vep_data=config['vep_data'],
        vep_synonyms=config['vep_synonyms'],
        gdc_fasta=config['genome_fasta'],
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcfVEP.txt"
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        "vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs --fa {params.gdc_fasta}"
    
rule vcf2maf:
    """General rule to convert the different vcf files into maf"""
    input:
        #NOTE: in vcf2maf.pl 1.6.18+ we can pass inhibit-vep which will
        #skip calling vep on the file
        vep="analysis/somatic/{run}/{run}_{caller}.{type}.vep.vcf",
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.maf"
    params:
        vep_index=config['vep_fasta'],
        vep_custom_enst= config['vep_custom_enst'],
        vep_assembly=config['vep_assembly'],
        vep_filter= config['vep_filter'],
        buffer_size=config['vcf2maf_bufferSize'],

        tumor= lambda wildcards: somatic_getTumor(wildcards),
        normal= lambda wildcards: somatic_getNormal(wildcards),
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcf2maf.txt"
    log:
        "analysis/logs/somatic/{run}/{run}.{caller}.{type}_vcf2maf.log.txt"
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        """vcf2maf.pl --input-vcf {input.vep} --output-maf {output} --custom-enst {params.vep_custom_enst} --ref-fasta {params.vep_index} --tumor-id {params.tumor} --normal-id {params.normal} --ncbi-build {params.vep_assembly} --filter-vcf {params.vep_filter} --buffer-size {params.buffer_size} --inhibit-vep 1 2> {log}"""


rule mutationSignature:
    """General rule to do mutation signature analysis using mutProfiler.py"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.{type}.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.{type}.pdf"
    params:
        index= lambda wildcards: os.path.abspath(config['genome_fasta']),
        #BUILD up our tcga_panel using a helper fn
        tcga_panel = build_tcga_param(),
        outname = lambda wildcards: "%sanalysis/somatic/%s/%s_%s.%s" % (config['remote_path'], wildcards.run, wildcards.run, wildcards.caller, wildcards.type),
        name = lambda wildcards: wildcards.run
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_mutationSignature.txt"
    group: "somatic"
    shell:
        "cidc_wes/cidc-vs/mutProfile.py -c {params.tcga_panel} -m {input} -r {params.index} -o {params.outname} -n {params.name}"

rule maf_exon_filter:
    """General rule to filter coding exon mutations"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_maf_exon_filter.txt"
    shell:
        "cidc_wes/modules/scripts/maf_exon_filter.py -m {input} -o {output}"

rule calculate_mutation:
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.mutationload.txt"
    params:
        size=config['effective_size'],
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.calculate_mutation.txt"
    shell:
        "cidc_wes/modules/scripts/mutation_load.py -v {input} -o {output} -s {params.size}"

rule extract_VAF_DEPTH:
    """Run Jingxins harmonization script to extract VAF and DEPTH stats"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.stats.txt"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.extract_VAF_DEPTH.txt"
    shell:
        "cidc_wes/modules/scripts/extract_vaf_depth.py -v {input} > {output}"

rule somatic_gzip_filtered_vcf:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz",
    group: "somatic"
    shell:
        "bgzip -c {input} > {output}"

rule somatic_tabix_filtered_vcf_gz:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz.tbi",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        "tabix -p vcf {input}"

rule somatic_getExonic_mutations:
    """Get the mutations that fall into the exonic regions"""
    input:
        vcf="analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz",
        tbi="analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz.tbi"
    params:
        exons=config['CDS_Bed_input']
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.somatic_getExonic_mutations.txt"
    shell:
        "bcftools view -R {params.exons} {input.vcf} | bcftools sort | bcftools view -Oz > {output}"

rule somatic_tabix_exonic_mutations:
    """Get the mutations that fall into the exonic regions"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz.tbi",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    shell:
        "tabix -p vcf {input}"

rule somatic_getTarget_mutations:
    """Get the mutations that fall into the exonic regions"""
    input:
        vcf="analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
        tbi="analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz.tbi"
    params:
        target= lambda wildcards: center_targets[wildcards.center]
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.exons.{center}.vcf.gz",
    group: "somatic"
    conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.{center}.somatic_getTarget_mutations.txt"
    shell:
        "bcftools view -R {params.target} {input} | bcftools sort | bcftools view -Oz > {output}"

rule summarize_somatic_mutations:
    """Use the filter.maf to generate summary statistics for SNPS, INS, DEL
    --used in the wes report"""
    input:
        expand("analysis/somatic/{run}/{run}_{{caller}}.filter.maf", run=sorted(config['runs']))
    output:
        cts = "analysis/somatic/somatic_mutation_summaries.{caller}.csv",
        annot = "analysis/somatic/somatic_functional_annot_summaries.{caller}.csv",
    params:
        files = lambda wildcards, input: " -m ".join(input)
    group: "somatic"
    benchmark:
        "benchmarks/somatic/summarize_somatic_mutations.{caller}.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_genStats.py -m {params.files} -o {output.cts} -a {output.annot}"

rule summarize_SNV_mutations:
    """Use the filter.vep.vcf to generate summary table for transition count
    table, e.g. outputs:
    Ref(rows)/Alt(cols),A,C,G,T
    A
    C
    G
    T
    --used in the wes report"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vep.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}_somatic_SNV_summaries.csv"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}_summarize_SNV_mutations.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_SNV_stats.py -v {input} -o {output}"

rule summarize_processINDELcircos:
    """Process the filter.maf file to generate a file suitable for circos"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.indel.circos.txt"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/summarize_processINDELcircos.{run}.{caller}.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_processINDEL.py -m {input} > {output}"

rule summarize_processSNPcircos:
    """Process the filter.maf file to generate a file suitable for circos"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.maf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.snp.circos.txt"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/summarize_processSNPcircos.{run}.{caller}.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_processSNP.py -m {input} > {output}"

rule somatic_json_filtered_maf:
    """json encode the filtered maf file; write out the filtered maf file as
    base64 string"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.maf"
    output:
        "analysis/report/json/somatic/{run}_{caller}.filtered_maf.json"
    params:
        run = lambda wildcards: wildcards.run
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}_{caller}.somatic_json_filtered_maf.txt"
    shell:
        "cidc_wes/modules/scripts/json_filtered_maf.py -r {params.run} -f {input} -o {output}"
