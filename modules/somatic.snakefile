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
    #center = config.get('cimac_center', 'broad') #Try to get center, default broad
    for run in config['runs']:
        ls.append("analysis/somatic/%s/%s_%s.output.vcf.gz" % (run,run, caller))
        #Optimized
        ls.append("analysis/somatic/%s/%s_%s.output.twist.vcf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.output.twist.maf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.output.twist.filtered.vcf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.output.twist.filtered.maf" % (run,run, caller))

        #Filtered
        ls.append("analysis/somatic/%s/%s_%s.filter.vcf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.vcf.gz" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.filter.maf" % (run,run, caller))
        #lego plot script, mutProfile.py generates a plot (.pdf) and a counts
        #file (.json)
        ls.append("analysis/somatic/%s/%s_%s.twist.pdf" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.twist.tri_mtrx.json" % (run,run, caller))
        #Sample summary report
        ls.append("analysis/somatic/%s/%s_%s_somatic_SNV_summaries.csv" % (run,run, caller))
        ls.append("analysis/somatic/%s/%s_%s.mutation_summaries.csv" % (run,run,caller))
        ls.append("analysis/somatic/%s/%s_%s.functional_annot_summaries.csv" % (run,run,caller))
        #Sample Json for cohort report
        ls.append("analysis/somatic/%s/%s_%s_onco_gene_list.tsv" % (run,run,caller))
    #json files
    ls.append("analysis/report/json/somatic/%s_%s.somatic.json" % (run, caller))
    return ls

def somatic_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    caller = config.get('somatic_caller', 'tnscope')
    ls = somatic_helper_targets(wildcards, caller)
        
    return ls

# Deprecated--no longer used!
# def somatic_output_files(wildcards):
#     """returns a list of filepaths generated by this module to store 
#     in the CIDC for a given sample 
#     """
#     ls = []
#     caller = config['somatic_caller']
#     #center = config.get('cimac_center', 'broad')
#     run = wildcards.run
#     ls.append("analysis/somatic/%s/%s_%s.filter.exons.center_targets.vcf.gz" % (run,run,caller))
#     ls.append("analysis/somatic/%s/%s_%s.filter.maf" % (run,run,caller))
#     ls.append("analysis/somatic/%s/%s_%s.filter.vcf" % (run,run,caller))
#     ls.append("analysis/somatic/%s/%s_%s.output.maf" % (run,run,caller))
#     ls.append("analysis/somatic/%s/%s_%s.output.vcf" % (run,run,caller))
    
#     return ls
# 
# def somatic_make_file_map_makeKeys():
#     """Makes the keys for the yaml_writer"""
#     center = config.get('cimac_center', 'broad')
#     tmp = " -k ".join(['%s_filter_exon' % center,'filter_maf','filter_vcf','output_maf','output_vcf']),

# rule somatic_make_file_map:
#     input:
#         somatic_output_files
#     output:
#         "analysis/somatic/{run}/{run}.somatic.output.yaml"
#     benchmark: "benchmarks/somatic/{run}/{run}.somatic_make_file_map.txt"
#     params:
#         run = lambda wildcards: wildcards.run,
#         #keys = " -k ".join(['broad_filter_exons','mda_filter_exons','mocha_filter_exon','filter_maf','filter_vcf','output_maf','output_vcf']),
#         kkeys = lambda wildcards: somatic_make_file_map_makeKeys(),
#         files = lambda wildcards, input: " -f ".join(input),
#     shell:
#         "cidc_wes/modules/scripts/yaml_writer.py -t runs -n {params.run} -k {params.kkeys} -f {params.files} > {output}"

rule somatic_all:
    input:
        somatic_targets
    benchmark: "benchmarks/somatic/somatic_all.txt"

###############################################################################
# PLEASE LOOK AT the different caller snakefiles,
# somatic_tnhaplotyper2, somatic_tnsnv, somatic_tnscope
# for caller specific rules!
###############################################################################

def somatic_twist_inputFn(wildcards):
    run = wildcards.run
    caller = wildcards.caller
    #default to this
    tmp = {'vcf': "analysis/somatic/%s/%s_%s.output.vcf.gz" % (run, run, caller),
           'tbi': "analysis/somatic/%s/%s_%s.output.vcf.gz.tbi" % (run, run, caller)}
    if config.get('tumor_only', False) and caller == 'tnscope':
        #For tumor only runs, when we use tnscope, we use preprocess.vcf.gz
        tmp = {'vcf': "analysis/somatic/%s/%s_preprocess.vcf.gz" % (run, run),
               'tbi': "analysis/somatic/%s/%s_preprocess.vcf.gz.tbi" % (run, run)}
    return tmp

rule somatic_twist:
    """Takes output.vcf.gz and intersects it with the twist capture regions"""
    input:
        unpack(somatic_twist_inputFn)
    output:
        "analysis/somatic/{run}/{run}_{caller}.output.twist.vcf",
    params:
        twist_regions = config['twist_regions']
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_somatic_twist.txt"
    shell:
        """bcftools view -R {params.twist_regions} {input.vcf} | bcftools sort | bcftools view -Ov > {output}"""

rule somatic_twist_filter:
    """Filters the PASS column of a vcf and removes evertying 
    except PASS and alt_allele_in_normal"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.twist.vcf",
    output:
        "analysis/somatic/{run}/{run}_{caller}.output.twist.filtered.vcf",
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_somatic_optimize.txt"    
    shell:
        """grep -v 'germline-risk\|low_t_alt_frac\|t_lod_fstar\|triallelic_site' {input} > {output}"""

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
	vcf_bin_path="%s/bin/" % config['vcf_root'],
    group: "somatic"
    #NOTE: b/c the rule uses a run instead of a shell, snkmk doesnt allow
    #conda env defs (next line)
    #conda: "../envs/somatic_vcftools.yml"
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}_filter_raw_vcf.txt"
    shell:
        """{params.vcf_bin_path}vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""
        #"""vcftools --gzvcf {input} --remove-filtered-all --recode --stdout > {output}"""


#LEN: Don't we need this rule??
# rule gunzip_vcf:
#     """General rule to gunzip the three types of vcf.gz files-
#     tnscope_, tnsnv, and tnhaplotyper"""
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.output.vcf.gz"
#     output:
#         #Should we make this a temp?
#         "analysis/somatic/{run}/{run}_{caller}.output.vcf"
#     benchmark:
#         "benchmarks/somatic/{run}/{run}.{caller}_gunzip_vcf.txt"
#     group: "somatic"
#     shell:
#         #NOTE: we want to keep the original .gz vcf file
#         "gunzip < {input} > {output}"

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
	vcf_bin_path="%s/bin/" % config['vcf_root'],
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcfVEP.txt"
    group: "somatic"
    shell:
        "{params.vcf_bin_path}vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs --fa {params.gdc_fasta} --format vcf"
        #"vep --i {input} --dir_cache={params.vep_data} --synonyms {params.vep_synonyms} --vcf -o {output} --offline --hgvs --fa {params.gdc_fasta} --format vcf"
    
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
	vcf_bin_path="%s/bin/" % config['vcf_root'],
        tumor= lambda wildcards: somatic_getTumor(wildcards),
        normal= lambda wildcards: somatic_getNormal(wildcards),
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.{type}_vcf2maf.txt"
    log:
        "analysis/logs/somatic/{run}/{run}.{caller}.{type}_vcf2maf.log.txt"
    group: "somatic"
    conda: "../envs/vcf.yml"
    shell:
        """{params.vcf_bin_path}vcf2maf.pl --input-vcf {input.vep} --output-maf {output} --custom-enst {params.vep_custom_enst} --ref-fasta {params.vep_index} --tumor-id {params.tumor} --normal-id {params.normal} --ncbi-build {params.vep_assembly} --filter-vcf {params.vep_filter} --buffer-size {params.buffer_size} --inhibit-vep 1 2> {log}"""
        #"""vcf2maf.pl --input-vcf {input.vep} --output-maf {output} --custom-enst {params.vep_custom_enst} --ref-fasta {params.vep_index} --tumor-id {params.tumor} --normal-id {params.normal} --ncbi-build {params.vep_assembly} --filter-vcf {params.vep_filter} --buffer-size {params.buffer_size} --inhibit-vep 1 2> {log}"""

rule mutationSignature:
    """General rule to do mutation signature analysis using mutProfiler.py"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.output.twist.maf"
    output:
        pdf="analysis/somatic/{run}/{run}_{caller}.twist.pdf",
        json="analysis/somatic/{run}/{run}_{caller}.twist.tri_mtrx.json",
    params:
        index= lambda wildcards: os.path.abspath(config['genome_fasta']),
        #BUILD up our tcga_panel using a helper fn
        tcga_panel = build_tcga_param(),
        outname = lambda wildcards: "%sanalysis/somatic/%s/%s_%s.twist" % (config['remote_path'], wildcards.run, wildcards.run, wildcards.caller),
        name = lambda wildcards: wildcards.run
    benchmark:
        "benchmarks/somatic/{run}/{run}.{caller}.mutationSignature.txt"
    group: "somatic"
    shell:
        "cidc_wes/cidc-vs/mutProfile.py -c {params.tcga_panel} -m {input} -r {params.index} -o {params.outname} -n {params.name} -j {output.json}"

#DEPRECATED b/c not used anymore
# rule maf_exon_filter:
#     """General rule to filter coding exon mutations"""
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.output.maf"
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
#     group: "somatic"
#     benchmark:
#         "benchmarks/somatic/{run}/{run}.{caller}_maf_exon_filter.txt"
#     shell:
#         "cidc_wes/modules/scripts/maf_exon_filter.py -m {input} -o {output}"

# rule calculate_mutation:
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.output.exon.maf"
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.mutationload.txt"
#     params:
#         size=config['effective_size'],
#     group: "somatic"
#     benchmark:
#         "benchmarks/somatic/{run}/{run}_{caller}.calculate_mutation.txt"
#     shell:
#         "cidc_wes/modules/scripts/mutation_load.py -v {input} -o {output} -s {params.size}"

# rule extract_VAF_DEPTH:
#     """Run Jingxins harmonization script to extract VAF and DEPTH stats"""
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.filter.stats.txt"
#     group: "somatic"
#     benchmark:
#         "benchmarks/somatic/{run}/{run}_{caller}.extract_VAF_DEPTH.txt"
#     shell:
#         "cidc_wes/modules/scripts/extract_vaf_depth.py -v {input} > {output}"

rule somatic_gzip_filtered_vcf:
    """Prepping the files filtered.vcf file for somatic_getExonic_mutations
    bgzip-ing and tabix"""
    input:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz",
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.gzip_filtered_vcf.txt"
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
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}.tabix_filtered_vcf_gz.txt"
    group: "somatic"
    #conda: "../envs/somatic_vcftools.yml"
    shell:
        "tabix -p vcf {input}"

# Deprecated b/c the targeted bed regions capture the exons we really should
# study
# rule somatic_getExonic_mutations:
#     """Get the mutations that fall into the exonic regions"""
#     input:
#         vcf="analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz",
#         tbi="analysis/somatic/{run}/{run}_{caller}.filter.vcf.gz.tbi"
#     params:
#         exons=config['CDS_Bed_input']
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
#     group: "somatic"
#     conda: "../envs/somatic_vcftools.yml"
#     benchmark:
#         "benchmarks/somatic/{run}/{run}_{caller}.getExonic_mutations.txt"
#     shell:
#         "bcftools view -R {params.exons} {input.vcf} | bcftools sort | bcftools view -Oz > {output}"

# rule somatic_tabix_exonic_mutations:
#     """Get the mutations that fall into the exonic regions"""
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz.tbi",
#     benchmark:
#         "benchmarks/somatic/{run}/{run}_{caller}.tabix_exonic_mutations.txt"
#     group: "somatic"
#     conda: "../envs/somatic_vcftools.yml"
#     shell:
#         "tabix -p vcf {input}"

#DEPRECATED b/c AFTER harmonization the targeted bed files are all integrated
#nto a sincle file - twist.broad.mdacc.mocha.liftover.hg38.sorted.merged.bed 
# rule somatic_getTarget_mutations:
#     """Get the mutations that fall into the exonic regions"""
#     input:
#         vcf="analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz",
#         tbi="analysis/somatic/{run}/{run}_{caller}.filter.exons.vcf.gz.tbi"
#     params:
#         target= lambda wildcards: center_targets[config.get('cimac_center', 'broad')]
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.filter.exons.center_targets.vcf.gz",
#     group: "somatic"
#     conda: "../envs/somatic_vcftools.yml"
#     benchmark:
#         "benchmarks/somatic/{run}/{run}_{caller}.getTarget_mutations.txt"
#     shell:
#         "bcftools view -R {params.target} {input} | bcftools sort | bcftools view -Oz > {output}"

rule summarize_somatic_mutations:
    """Use the filter.maf to generate summary statistics for SNPS, INS, DEL
    --used in the wes report"""
    input:
        #expand("analysis/somatic/{run}/{run}_{{caller}}.filter.maf", run=sorted(config['runs'])) #NO LONGER expecting multiple runs!
        maf="analysis/somatic/{run}/{run}_{caller}.output.twist.maf"
    output:
        cts = "analysis/somatic/{run}/{run}_{caller}.mutation_summaries.csv",
        annot = "analysis/somatic/{run}/{run}_{caller}.functional_annot_summaries.csv",
    params:
        #files = lambda wildcards, input: " -m ".join(input),
        targets = center_targets[config.get('cimac_center', 'broad')],
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_summarize_somatic_mutations.{caller}.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_genStats.py -m {input.maf} -t {params.targets} -o {output.cts} -a {output.annot}"

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
        "analysis/somatic/{run}/{run}_{caller}.output.twist.vep.vcf"
    output:
        "analysis/somatic/{run}/{run}_{caller}_somatic_SNV_summaries.csv"
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}_summarize_SNV_mutations.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_SNV_stats.py -v {input} -o {output}"

# DEPRECATED b/c no used?
# rule summarize_processINDELcircos:
#     """Process the filter.maf file to generate a file suitable for circos"""
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.filter.maf"
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.indel.circos.txt"
#     group: "somatic"
#     benchmark:
#         "benchmarks/somatic/summarize_processINDELcircos.{run}.{caller}.txt"
#     shell:
#         "cidc_wes/modules/scripts/somatic_processINDEL.py -m {input} > {output}"
#
# rule summarize_processSNPcircos:
#     """Process the filter.maf file to generate a file suitable for circos"""
#     input:
#         "analysis/somatic/{run}/{run}_{caller}.filter.maf"
#     output:
#         "analysis/somatic/{run}/{run}_{caller}.snp.circos.txt"
#     group: "somatic"
#     benchmark:
#         "benchmarks/somatic/summarize_processSNPcircos.{run}.{caller}.txt"
#     shell:
#         "cidc_wes/modules/scripts/somatic_processSNP.py -m {input} > {output}"

###############################################################################
# FOR cohort report
###############################################################################
rule somatic_json:
    """include the twist maf file (base64 encoded), trinucleotide matrix, 
    and TMB via summaries table"""
    input:
        maf="analysis/somatic/{run}/{run}_{caller}.output.twist.maf",
        tri_mtrx="analysis/somatic/{run}/{run}_{caller}.twist.tri_mtrx.json",
        summary="analysis/somatic/{run}/{run}_{caller}.mutation_summaries.csv",
        oncoGeneList="analysis/somatic/{run}/{run}_{caller}_onco_gene_list.tsv",
    output:
        "analysis/report/json/somatic/{run}_{caller}.somatic.json"
    params:
        run = lambda wildcards: wildcards.run
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}_{caller}.somatic_json.txt"
    shell:
        "cidc_wes/modules/scripts/json_somatic.py -r {params.run} -m {input.maf} -j {input.tri_mtrx} -s {input.summary} -l {input.oncoGeneList} -o {output}"

rule somatic_get_top_oncogenes:
    """Use the cancerGeneList.tsv and check to see which ones are represented
    in the twist maf file"""
    input:
        maf="analysis/somatic/{run}/{run}_{caller}.output.twist.maf",
        cancerGeneList = "cidc_wes/static/oncoKB/cancerGeneList.tsv",
    output:
        "analysis/somatic/{run}/{run}_{caller}_onco_gene_list.tsv",
    #params:
        #targets = center_targets[config.get('cimac_center', 'broad')],
    group: "somatic"
    benchmark:
        "benchmarks/somatic/{run}/{run}_{caller}_somatic_get_top_oncogenes.txt"
    shell:
        "cidc_wes/modules/scripts/somatic_getTopOncoGenes.py -m {input.maf} -l {input.cancerGeneList} -o {output}"
