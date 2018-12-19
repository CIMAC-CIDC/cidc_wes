import "../cidc-pipelines/wdl/wes/alignment.wdl" as alignment_wdl
import "../cidc-pipelines/wdl/wes/mutect2.wdl" as mutect2_wdl
import "../cidc-pipelines/wdl/wes/vep_offline.wdl" as vep_offline_wdl
import "../cidc-pipelines/wdl/preprocessing/fastqc.wdl" as fastqc_wdl
import "../cidc-pipelines/wdl/wes/selectvariants.wdl" as selectvariants_wdl
import "../cidc-pipelines/wdl/wes/funcotator.wdl" as funcotator_wdl

workflow run_wes {
    File tumor_fastq1
    File tumor_fastq2
    File normal_fastq1
    File normal_fastq2
    String prefix
    
    # BWA requirements
    Array[File] bwa_index

    Int bwa_memory
    Int bwa_disk_space
    Int bwa_num_cpu
    Int bwa_num_preempt
    Int bwa_boot_disk_gb

    # GATK requirements
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dictionary
    File KnownSites_dbsnp
    File KnownSites_dbsnp_index
    File KnownSites_indels
    File KnownSites_indels_index

    call fastqc_wdl.run_fastqc as tumor_fastqc1 {
        input:
            fastq = tumor_fastq1,
            prefix = prefix+'.T1',
            type = 'left'
    }
    call fastqc_wdl.run_fastqc as tumor_fastqc2 {
        input:
            fastq = tumor_fastq2,
            prefix = prefix+'.T2',
            type = 'right'
    }
    call fastqc_wdl.run_fastqc as normal_fastqc1 {
        input:
            fastq = normal_fastq1,
            prefix = prefix+'.N1',
            type = 'left'
    }
    call fastqc_wdl.run_fastqc as normal_fastqc2 {
        input:
            fastq = normal_fastq2,
            prefix = prefix+'.N2',
            type = 'right'
    }

    call alignment_wdl.run_alignment as tumor_align {
        input:
            fastq1 = tumor_fastq1,
            fastq2 = tumor_fastq2,
            prefix = prefix+'.tumor',
            groupname = "tumor",
            bwa_index = bwa_index,
            genome_fasta = genome_fasta,
            genome_fasta_index = genome_fasta_index,
            genome_fasta_dictionary = genome_fasta_dictionary,
            KnownSites_dbsnp = KnownSites_dbsnp,
            KnownSites_dbsnp_index = KnownSites_dbsnp_index,
            KnownSites_indels = KnownSites_indels,
            KnownSites_indels_index = KnownSites_indels_index,
            bwa_memory = bwa_memory,
            bwa_disk_space = bwa_disk_space,
            bwa_num_preempt = bwa_num_preempt,
            bwa_num_cpu = bwa_num_cpu,
            bwa_boot_disk_gb = bwa_boot_disk_gb
    }
    call alignment_wdl.run_alignment as normal_align {
        input:
            fastq1 = normal_fastq1,
            fastq2 = normal_fastq2,
            prefix = prefix+'.normal',
            groupname = "normal",
            bwa_index = bwa_index,
            genome_fasta = genome_fasta,
            genome_fasta_index = genome_fasta_index,
            genome_fasta_dictionary = genome_fasta_dictionary,
            KnownSites_dbsnp = KnownSites_dbsnp,
            KnownSites_dbsnp_index = KnownSites_dbsnp_index,
            KnownSites_indels = KnownSites_indels,
            KnownSites_indels_index = KnownSites_indels_index,
            bwa_memory = bwa_memory,
            bwa_disk_space = bwa_disk_space,
            bwa_num_preempt = bwa_num_preempt,
            bwa_num_cpu = bwa_num_cpu,
            bwa_boot_disk_gb = bwa_boot_disk_gb
    }
    call mutect2_wdl.run_mutect2 as run_mutect2 {
        input:
            tumor_bam = tumor_align.recalibrated_bam,
            tumor_bam_index = tumor_align.recalibrated_bam_index,
            normal_bam = normal_align.recalibrated_bam,
            normal_bam_index = normal_align.recalibrated_bam_index,
            ref_fasta = genome_fasta,
            ref_fasta_index = genome_fasta_index,
            ref_dict = genome_fasta_dictionary,
            prefix = prefix,
            read_groups_set = true
    }
    call selectvariants_wdl.selectvariants as selectvariants {
        input:
            input_vcf = run_mutect2.filtered_vcf,
            input_vcf_index = run_mutect2.filtered_vcf_index
    }
    call vep_offline_wdl.run_vep as run_vep {
        input:
            inputFile = selectvariants.output_vcf,
            prefix = prefix,
            ref_fasta = genome_fasta
    }
    call funcotator_wdl.run_funcotator as run_funcotator {
        input:
            input_vcf = run_mutect2.filtered_vcf,
            input_vcf_index = run_mutect2.filtered_vcf_index,
            prefix = prefix,
            genome_fasta = genome_fasta,
            genome_fasta_index = genome_fasta_index,
            genome_fasta_dict = genome_fasta_dictionary
    }
    output {
        #FastQC outputs
        File normal_1_fastqc_html = normal_fastqc1.html
        File normal_1_fastqc_zip = normal_fastqc1.zip
        Map[String, String] normal_1_fastqc_summary = normal_fastqc1.summary
        File normal_2_fastqc_html = normal_fastqc2.html
        File normal_2_fastqc_zip = normal_fastqc2.zip
        Map[String, String] normal_2_fastqc_summary = normal_fastqc2.summary
        File tumor_1_fastqc_html = tumor_fastqc1.html
        File tumor_1_fastqc_zip = tumor_fastqc1.zip
        Map[String, String] tumor_1_fastqc_summary = tumor_fastqc1.summary
        File tumor_2_fastqc_html = tumor_fastqc2.html
        File tumor_2_fastqc_zip = tumor_fastqc2.zip
        Map[String, String] tumor_2_fastqc_summary = tumor_fastqc2.summary

        # Alignment outputs
        File normal_align_recalibrated_bam = normal_align.recalibrated_bam
        File normal_align_recalibrated_bam_index = normal_align.recalibrated_bam_index
        File tumor_align_recalibrated_bam = tumor_align.recalibrated_bam
        File tumor_align_recalibrated_bam_index = tumor_align.recalibrated_bam_index

        # Mutect2 outputs
        File run_mutect2_unfiltered_vcf = run_mutect2.unfiltered_vcf
        File run_mutect2_unfiltered_vcf_index = run_mutect2.unfiltered_vcf_index
        File run_mutect2_filtered_vcf = run_mutect2.filtered_vcf
        File run_mutect2_filtered_vcf_index = run_mutect2.filtered_vcf_index

        # Variants selected for VEP
        File selectedvariants_for_vep_vcf = selectvariants.output_vcf
        File selectedvariants_for_vep_vcf_index = selectvariants.output_vcf_index

        # VEP outputs
        File vep_annotated_vcf = run_vep.output_vep_vcf 
        File vep_annotated_maf = run_vep.output_vep_maf

        # Funcotator outputs
        File funcotator_annotated_vcf = run_funcotator.funcotator_output_vcf
        File funcotator_annotated_maf = run_funcotator.funcotator_output_maf
    }
    meta {
        maintainer: "Jason L Weirather"
        citation1: "Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN]"
        citation2: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
        citation3: "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002."
    }
}

