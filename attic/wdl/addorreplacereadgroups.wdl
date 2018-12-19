task addorreplacereadgroups {
    File input_bam
    String? prefix = "addrepgroup"
    String? RGLB = "lib1"
    String? RGPL = "illumina"
    String? RGPU = "1"
    String? RGSM = "1"
    String? RGID = "1"


    # Required options
    Int? memory = 16
    Int? disk_space = 100
    Int? num_cpu = 2
    Int? num_preempt = 0
    Int? boot_disk_gb = 50

    String? docker = "broadinstitute/gatk:4.0.3.0"

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
       gatk --java-options "-Xmx4g"  AddOrReplaceReadGroups \
           --INPUT ${input_bam} \
           --OUTPUT ${prefix}.readgroups.bam \
           --RGLB ${RGLB} \
           --RGPL ${RGPL} \
           --RGPU ${RGPU} \
           --RGSM ${RGSM} \
           --RGID ${RGID} \
           --VALIDATION_STRINGENCY LENIENT 2> /dev/null
    >>>
    output {
       File output_bam = "${prefix}.readgroups.bam"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
    }
}

workflow run_addorreplacereadgroups {
    call addorreplacereadgroups {
    }
}
