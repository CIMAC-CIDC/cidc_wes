task variant_effect_predictor {
    File inputFile
    String outputFileName
    command {
     vep \
        --i ${inputFile}  \
        --o ${outputFileName} --database

        }

    output {
        File annotatedFile = "${outputFileName}"
    }

    runtime {
        docker: "ensemblorg/ensembl-vep"
    }
}

workflow vep {
    call variant_effect_predictor
}
