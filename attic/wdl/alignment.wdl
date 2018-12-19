import "../cidc-pipelines/wdl/wes/bwa.wdl" as bwa_wdl
import "../cidc-pipelines/wdl/wes/gatk_preprocessing.wdl" as gatk_preprocessing_wdl
import "../cidc-pipelines/wdl/utilities/idxstats.wdl" as idx

workflow run_alignment {
    File fastq1
    File fastq2
    String prefix
    String groupname
    
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

    call bwa_wdl.bwamem as bwamem {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            index = bwa_index,
            prefix = prefix,
            disk_space = bwa_disk_space,
            memory = bwa_memory,
            num_cpu = bwa_num_cpu,
            sort_threads = 4,
            complevel = 1
    }
    call idx.idxstats as idxstats {
       input: bam=bwamem.bam,
              bam_index=bwamem.bam_index
    }
    call gatk_preprocessing_wdl.run_gatk_preprocessing as run_gatk_preprocessing {
        input:
            input_bam = bwamem.bam,
            input_bam_index = bwamem.bam_index,            
            prefix = prefix,
            groupname = groupname,
            genome_fasta = genome_fasta,
            genome_fasta_index = genome_fasta_index, 
            genome_fasta_dictionary = genome_fasta_dictionary,
            KnownSites_dbsnp = KnownSites_dbsnp,
            KnownSites_dbsnp_index = KnownSites_dbsnp_index,
            KnownSites_indels = KnownSites_indels,
            KnownSites_indels_index = KnownSites_indels_index
    }
    output {
        Array[Array[String]] raw_idxstats_table = idxstats.idxstats_table
        Float raw_fraction_aligned_reads = idxstats.fraction_aligned_reads
        Int raw_total_reads = idxstats.total_reads
        Int raw_aligned_reads = idxstats.aligned_reads
        Int raw_unaligned_reads = idxstats.unaligned_reads
        File raw_bam = bwamem.bam
        File raw_bam_index = bwamem.bam_index

        Array[Array[String]] recalibration_idxstats_table = run_gatk_preprocessing.idxstats_table
        Float recalibration_fraction_aligned_reads = run_gatk_preprocessing.fraction_aligned_reads
        Int recalibration_total_reads = run_gatk_preprocessing.total_reads
        Int recalibration_aligned_reads = run_gatk_preprocessing.aligned_reads
        Int recalibration_unaligned_reads = run_gatk_preprocessing.unaligned_reads
        File recalibrated_bam = run_gatk_preprocessing.recalibrated_bam
        File recalibrated_bam_index = run_gatk_preprocessing.recalibrated_bam_index
    }
    meta {
        maintainer: "Jason L Weirather"
        citation1: "Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN]"
        citation2: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
        citation3: "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002."
    }
}

