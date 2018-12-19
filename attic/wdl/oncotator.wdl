# Get Mutect2 hg38 VCF file
# Liftover to hg19  VCF file
# Run Oncotator in hg19

workflow Oncotator {
    File input_vcf
    File hg38_hg19_chain
    File hg19_fasta
    File? hg19_dict

    File default_config_file
    String? prefix="output"
    File? onco_ds_tar_gz
    String? oncotator_exe
    String? sequencing_center
    String? sequence_source
    String? control_id

    # Required options
    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    call liftover_vcf {
       input: input_vcf=input_vcf,
              liftover_chain=hg38_hg19_chain,
              reference_fasta=hg19_fasta,
              reference_dict=hg19_dict,
              prefix=prefix,
              memory=memory,
              disk_space=disk_space,
              num_cpu=num_cpu,
              num_preempt=num_preempt,
              boot_disk_gb=boot_disk_gb
    }

    call oncotate_m2 {
       input: m2_vcf=liftover_vcf.valid_vcf,
              case_id=prefix,
              onco_ds_tar_gz=onco_ds_tar_gz,
              oncotator_exe=oncotator_exe,
              sequencing_center=sequencing_center,
              sequence_source=sequence_source,
              default_config_file=default_config_file,
              control_id=control_id
    }

    output {
       File oncotated_maf=oncotate_m2.oncotated_m2_maf
    }
}

task liftover_vcf {
   File input_vcf
   File liftover_chain
   File reference_fasta
   File? reference_dict
   String? prefix = "output"

   # Required options
   Int? memory = 8
   Int? disk_space = 20
   Int? num_cpu = 4
   Int? num_preempt = 0
   Int? boot_disk_gb = 10
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
   command {
        if [[ "${input_vcf}" == *.gz ]]; then
            echo "Decompress input vcf files"
            gunzip -c ${input_vcf} > decompressed.vcf
            vcf_file=decompressed.vcf
        else
            vcf_file=${input_vcf}
        fi

        if [[ $(ls ${reference_dict}) == '' ]]
        then
        java -Xmx8g -jar /gatk/gatk.jar CreateSequenceDictionary \
             -REFERENCE=${reference_fasta}
        fi

        java -Xmx8g -jar /gatk/gatk.jar LiftoverVcf \
              -I=$vcf_file \
              -O=${prefix}_valid.vcf \
          -CHAIN=${liftover_chain} \
         -REJECT=${prefix}_reject.vcf \
              -R=${reference_fasta}
   }
   output {
        File valid_vcf = "${prefix}_valid.vcf"
        File reject_vcf = "${prefix}_reject.vcf"
    }
}

task oncotate_m2 {
    File m2_vcf
    File default_config_file
    String? case_id="output"
    File? onco_ds_tar_gz
    String? oncotator_exe
    String? sequencing_center
    String? sequence_source
    String? control_id

    # Runtime parameters
    Int? mem = 10
    Int? preemptible_attempts = 2
    Int? disk_space_gb = 100
    String? docker = "broadinstitute/oncotator:1.9.6.1"

    command <<<

          set -e

          BASE=tempDir
          if [[ -d /mnt/ramdisk ]]; then
              BASE=/mnt/ramdisk/tempDir
          fi

          if [[ "${onco_ds_tar_gz}" == *.tar.gz ]]; then
              echo "Using given tar file: ${onco_ds_tar_gz}"
              mkdir $BASE
              tar zxf ${onco_ds_tar_gz} -C $BASE/ --strip-components 1
          else
              echo "Downloading and installing oncotator datasources from Broad FTP site..."
              # Download and untar the db-dir
              wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz
              tar zxf oncotator_v1_ds_April052016.tar.gz -C $BASE/  --strip-components 1
          fi

          chmod 777 -R $BASE
          ln -s $BASE onco_dbdir

        ${default="/root/oncotator_venv/bin/oncotator" oncotator_exe} --db-dir onco_dbdir/ -c $HOME/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt  \
            -v ${m2_vcf} ${case_id}.maf.annotated hg19 -i VCF -o TCGAMAF --skip-no-alt --infer-onps --collapse-number-annotations --log_name oncotator.log \
            -a Center:${default="Unknown" sequencing_center} \
            -a source:${default="Unknown" sequence_source} \
            -a normal_barcode:${control_id} \
            -a tumor_barcode:${case_id} \
            ${"--default_config " + default_config_file}
    >>>

    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: select_first([mem, 3]) + " GB"
        slurm_memory: select_first([mem, 3]) + "G"
        bootDiskSizeGb: 12
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
        zones: "us-east1-b us-east1-c us-east1-d"
        mount_tmpfs: "/mnt/ramdisk"
    }

    output {
        File oncotated_m2_maf="${case_id}.maf.annotated"
    }
}
