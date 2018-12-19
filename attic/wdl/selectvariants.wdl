task selectvariants {
    File input_vcf
    File input_vcf_index
    String? prefix = "variants"

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
       mkdir mytemp
       gatk  --java-options "-Xmx8g" SelectVariants \
           -V ${input_vcf} \
           -O ${prefix}.selected.vcf \
           --exclude-filtered
    >>>
    output {
       File output_vcf = "${prefix}.selected.vcf"
       File output_vcf_index = "${prefix}.selected.vcf.idx"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
    }
}

workflow run_selectvariants {
    call selectvariants {
    }
}
