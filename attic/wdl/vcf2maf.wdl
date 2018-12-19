task vcf2maf {
    File input_vcf
    File filter_vcf
    File ref_fasta
    String? tumor_label = "TUMOR"
    String? normal_label = "NORMAL"

    # Required options
    Int? memory = 16
    Int? disk_space = 100
    Int? num_cpu = 2
    Int? num_preempt = 0
    Int? boot_disk_gb = 20

    String? docker = "vacation/vcf2maf:1.6.16"
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
    command {
       cp "${input_vcf}" temp.vcf
       cp temp.vcf temp.vep.vcf
       perl /opt/vcf2maf/vcf2maf.pl --input-vcf temp.vcf  --output-maf "${tumor_label}.vep.maf" --custom-enst /opt/vcf2maf/data/isoform_overrides_uniprot --ref-fasta "${ref_fasta}" --tumor-id ${tumor_label} --normal-id ${normal_label} --ncbi-build GRCh38 --filter-vcf "${filter_vcf}"
    }
    output {
       File output_maf = "${tumor_label}.vep.maf"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "MSK vcf2maf https://github.com/mskcc/vcf2maf"
    }
}

workflow run_vcf2maf {
    File input_vcf
    File filter_vcf
    File ref_fasta
    String? tumor_label = "TUMOR"
    String? normal_label = "NORMAL"

    call vcf2maf {
        input:
            input_vcf = input_vcf,
            filter_vcf = filter_vcf,
            ref_fasta = ref_fasta,
            tumor_label = tumor_label,
            normal_label = normal_label
    }
}

