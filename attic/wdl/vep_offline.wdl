import "../cidc-pipelines/wdl/wes/vcf2maf.wdl" as vcf2maf_wdl

task vep {
    File inputFile
    File vep_cacheDir
    String? outputFileName = "output"
    Boolean? use_synonyms = true
    Boolean? vcf = true
    Int? disk_space = 200
    Int? boot_disk_gb = 40
    Int? memory = 32
    Int? num_cpu = 2
    Int? num_preempt = 0
    String? docker = "ensemblorg/ensembl-vep:release_91.3"
    runtime {
        docker: "${docker}"
        disks: "local-disk ${disk_space} SSD"
        slurm_docker: "${docker}"
        memory: "${memory} GB"
        slurm_memory: "${memory}G"
        cpu: "${num_cpu}"
        preemptible: "${num_preempt}"
        bootDiskSizeGb: "${boot_disk_gb}"
        zones: "us-east1-b us-east1-c us-east1-d"
        mount_tmpfs: "/mnt/ramdisk"
    }

    command {
        BASE=$(pwd)/tempDir
        if [ -d /mnt/ramdisk ]
          then
            BASE=/mnt/ramdisk/tempDir
        fi
        mkdir $BASE
        tar -xf ${vep_cacheDir} -C $BASE/.

        vep \
        --i ${inputFile}  \
        --dir_cache=$BASE/vep_cache \
        ${true=" --synonyms $BASE/vep_cache/homo_sapiens/91_GRCh38/chr_synonyms.txt" false="" use_synonyms} \
        ${true=" --vcf " false="" vcf} \
        --o "${outputFileName}" --offline --hgvs
        if ! [[ -d /mnt/ramdisk ]];
           then
              rm -r $(pwd)/tempDir
        fi
    }

    output {
        File annotatedFile = "${outputFileName}"
    }
}

workflow run_vep {
    File inputFile
    File vep_cacheDir
    String prefix
    File filter_vcf
    File ref_fasta
    String? vep_vcf_filename = "${prefix}.vep.vcf"
    Boolean? use_synonyms = true
    call vep {
        input:
            inputFile = inputFile,
            vep_cacheDir = vep_cacheDir,
            outputFileName = vep_vcf_filename,
            use_synonyms = use_synonyms
    }
    call vcf2maf_wdl.vcf2maf as vcf2maf {
        input:
            input_vcf = vep.annotatedFile,
            filter_vcf = filter_vcf,
            ref_fasta = ref_fasta,
            tumor_label = prefix
    }
    output {
        File output_vep_vcf = vep.annotatedFile
        File output_vep_maf = vcf2maf.output_maf
    }
}
