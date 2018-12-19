task baserecalibrator {
    File input_bam
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dictionary
    File KnownSites_dbsnp
    File KnownSites_dbsnp_index
    File KnownSites_indels
    File KnownSites_indels_index
    String? prefix = "baserecal"

    # Required options
    Int? memory = 16
    Int? disk_space = 100
    Int? num_cpu = 4
    Int? num_preempt = 0
    Int? boot_disk_gb = 20

    String? docker = "broadinstitute/gatk:4.0.1.2"

    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: "${memory} GB"
        slurm_memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${num_preempt}"
        bootDiskSizeGb: "${boot_disk_gb}"
        zones: "us-east1-b us-east1-c us-east1-d"
    }
    parameter_meta {
       input_bam: "bam in the proper format for processing by picard markduplicates suchas bwa mem -M"
    }
    command <<<

       gatk --java-options "-Xmx4g" BaseRecalibrator \
           -I ${input_bam} \
           -O recal_data.table \
           --reference  ${genome_fasta} \
           --known-sites ${KnownSites_dbsnp}  \
           --known-sites ${KnownSites_indels}

       gatk --java-options "-Xmx4g" ApplyBQSR \
           --bqsr-recal-file recal_data.table \
           -R ${genome_fasta} \
           -I ${input_bam} \
           -O ${prefix}.recalibrate.bam
    >>>
    output {
       File recalibration_data_table = "recal_data.table"
       File output_bam = "${prefix}.recalibrate.bam"
       File output_bai = "${prefix}.recalibrate.bai"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
    }
}

workflow run_baserecalibrator {
    call baserecalibrator {
    }
}
