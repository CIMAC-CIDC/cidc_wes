task variant_effect_predictor {
    File inputFile
    File? vep_cacheDir
    String outputFileName
    Int? use_online = 0
    Int disk_space
    Int? boot_disk_gb = 25
    String? docker = "ensemblorg/ensembl-vep"

    runtime {
        docker: "${docker}"
        disks: "local-disk ${disk_space} HDD"
        slurm_docker: "${docker}"
        bootDiskSizeGb: "${boot_disk_gb}"
        zones: "us-east1-b us-east1-c us-east1-d"
    }

    command {
    if [ ${use_online}==1 ]
    then
     vep \
        --i ${inputFile}  \
        --o ${outputFileName} --database

    else
        mkdir tempDir

        tar -xf ${vep_cacheDir} -C tempDir/.

        vep \
        --i ${inputFile}  \
        --dir_cache=tempDir/vep_cache \
        --o ${outputFileName} --offline --hgvs

        rm -rf tempDir
    fi

    }

    output {
        File annotatedFile = "${outputFileName}"
    }


}

workflow vep {
    call variant_effect_predictor
}
