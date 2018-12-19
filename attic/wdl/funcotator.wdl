task funcotator {
    File input_vcf
    File input_vcf_index
    File funcotator_data_tgz
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dict
    String reference_version # hg19 or hg38
    String? prefix = "funcotator"
    String output_format

    # Required options
    Int? memory = 32
    Int? disk_space = 100
    Int? num_cpu = 4
    Int? num_preempt = 0
    Int? boot_disk_gb = 20

    String? docker = "vacation/gatk:4.0.3.0-54295ab"

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
    }
    command <<<
       tar -xzf ${funcotator_data_tgz}
       mv funcotator_dataSources* funcotator_dataSources
       OUTFORMAT=MAF
       if [ "${output_format}" == "vcf" ]
       then
           OUTFORMAT=VCF
       fi
       gatk --java-options "-Xmx16g" Funcotator \
           --data-sources-path $(pwd)/funcotator_dataSources/ \
           --output-file-format $OUTFORMAT \
           --output "${prefix}.funcotator.${output_format}" \
           --ref-version ${reference_version} \
           -R ${genome_fasta}  \
           -V "${input_vcf}"
       #rm -r funcotator_dataSources/
    >>>
    output {
       File funcotator_output = "${prefix}.funcotator.${output_format}"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
    }
}

workflow run_funcotator {
    File input_vcf
    File input_vcf_index
    File funcotator_data_tgz
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dict
    String reference_version # hg19 or hg38
    String? prefix = "funcotator"
    call funcotator as vcf_funcotator {
        input: 
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            funcotator_data_tgz = funcotator_data_tgz,
            genome_fasta = genome_fasta,
            genome_fasta_index = genome_fasta_index,
            genome_fasta_dict = genome_fasta_dict,
            reference_version = reference_version,
            output_format="vcf",
            prefix = prefix
    }
    call funcotator as maf_funcotator {
        input: 
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            funcotator_data_tgz = funcotator_data_tgz,
            genome_fasta = genome_fasta,
            genome_fasta_index = genome_fasta_index,
            genome_fasta_dict = genome_fasta_dict,
            reference_version = reference_version,
            output_format="maf",
            prefix = prefix
    }
    output {
       File funcotator_output_vcf = vcf_funcotator.funcotator_output
       File funcotator_output_maf = maf_funcotator.funcotator_output
    }
}
