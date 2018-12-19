import "../cidc-pipelines/wdl/utilities/split_bam_sequences.wdl" as split_bam_sequences_wdl
import "../cidc-pipelines/wdl/wes/markduplicates.wdl" as markduplicates_wdl
import "../cidc-pipelines/wdl/wes/addorreplacereadgroups.wdl" as addorreplacereadgroups_wdl
import "../cidc-pipelines/wdl/wes/baserecalibrator.wdl" as baserecalibrator_wdl
import "../cidc-pipelines/wdl/utilities/merge_bams.wdl" as merge_bams_wdl
import "../cidc-pipelines/wdl/utilities/idxstats.wdl" as idxstats_wdl

workflow run_gatk_preprocessing {
    File input_bam
    File input_bam_index
    String? prefix = "gatk_preprocessed"
    String? groupname = "1"

    # base recalibrator parameters
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dictionary
    File KnownSites_dbsnp
    File KnownSites_dbsnp_index
    File KnownSites_indels
    File KnownSites_indels_index

    call split_bam_sequences_wdl.run_split_bam_sequences as run_split_bam_sequences {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index
    }
    ## Run the operations for everything except the unaligned reads
    scatter (bam_map in run_split_bam_sequences.output_split_bam_sequences) {
        call markduplicates_wdl.markduplicates as markduplicates {
            input:
                prefix = basename(bam_map["bam"],".bam"),
                input_bam = bam_map["bam"]
        }
        call addorreplacereadgroups_wdl.addorreplacereadgroups as addorreplacereadgroups {
            input:
                prefix = basename(bam_map["bam"],".bam"),
                input_bam = markduplicates.output_bam,
                RGSM = groupname,
                RGID = groupname,
                RGLB = groupname,
                RGPU = groupname,
                RGPL = groupname
        }
        call baserecalibrator_wdl.baserecalibrator as baserecalibrator {
            input:
                prefix = basename(bam_map["bam"],".bam"),
                input_bam = addorreplacereadgroups.output_bam,
                genome_fasta = genome_fasta,
                genome_fasta_index = genome_fasta_index,
                genome_fasta_dictionary = genome_fasta_dictionary,
                KnownSites_dbsnp = KnownSites_dbsnp,
                KnownSites_dbsnp_index = KnownSites_dbsnp_index,
                KnownSites_indels = KnownSites_indels,
                KnownSites_indels_index =KnownSites_indels_index
        }
    }
    ## The unaligned reads get their own special run through
    ##   skip base recalibrator on unaligned reads
    call markduplicates_wdl.markduplicates as markduplicatesUN {
            input:
                prefix = "unaligned",
                input_bam = run_split_bam_sequences.output_unaligend_bam_sequences["bam"]
    }
    call addorreplacereadgroups_wdl.addorreplacereadgroups as addorreplacereadgroupsUN {
            input:
                prefix = "unaligned",
                input_bam = markduplicatesUN.output_bam
    }
    call merge_bams_wdl.merge_bams as merge_bams {
        input:
            prefix = prefix+'.recalibrated',
            extra_bams = [addorreplacereadgroupsUN.output_bam],
            input_bams = baserecalibrator.output_bam
    }
    call idxstats_wdl.idxstats as recalibrated_idxstats {
        input:
            bam = merge_bams.output_bam,
            bam_index = merge_bams.output_bam_index
    }
    output {
        File recalibrated_bam = merge_bams.output_bam
        File recalibrated_bam_index = merge_bams.output_bam_index
        Array[Array[String]] idxstats_table = recalibrated_idxstats.idxstats_table
        Float fraction_aligned_reads = recalibrated_idxstats.fraction_aligned_reads
        Int total_reads = recalibrated_idxstats.total_reads
        Int aligned_reads = recalibrated_idxstats.aligned_reads
        Int unaligned_reads = recalibrated_idxstats.unaligned_reads
    }
    meta {
        maintainer: "Jason L Weirather"
        citation1: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
        citation2: "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002."
    }
}

