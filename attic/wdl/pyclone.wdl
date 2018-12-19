task pyclone {
    String prefix
    Array[File] pyclone_tsv
    String? docker="gcr.io/cidc-biofx/pyclone:0.13.1"
    Int? memory=8

    runtime {
        docker: "${docker}"
        slurm_docker: "${docker}"
        memory: "${memory} GB"
        slurm_memory: "${memory}G"
    }
    command {

        work=$PWD
        PyClone run_analysis_pipeline --in_files ${sep=' ' pyclone_tsv} --working_dir $work/${prefix}_pyclone_analysis
        tar -zcf $work/${prefix}_pyclone_analysis.tar.gz $work/${prefix}_pyclone_analysis

    }
    output {
        File pyclone_tar = "${prefix}_pyclone_analysis.tar.gz"
    }
    meta {
        maintainer: "Xihao Hu"
        citation: "Roth, A., Khattra, J., Yap, D., Wan, A., Laks, E., Biele, J., Ha, G., Aparicio, S., Bouchard-Côté, A. and Shah, S.P., 2014. PyClone: statistical inference of clonal population structure in cancer. Nature methods, 11(4), p.396."
    }
}

workflow PyClone {
    String prefix
    Array[File] pyclone_tsv
    Int? memory=8

    call pyclone {
       input: prefix=prefix,
              pyclone_tsv=pyclone_tsv,
              memory=memory
    }
}

