task seq_to_gistic {
    Array[File] seq_files

    command {
        mkdir inputs
        cp ${sep=' ' seq_files} inputs
        Rscript /sequenzaToGistic.r -i inputs -o .
    }

    output {
        File seq_file='all_segments.txt' 
        File markers_file='markers_gistic.txt'
    }

    runtime {
        docker: "gcr.io/cidc-biofx/gistic_r:v1"
        slurm_docker: "gcr.io/cidc-biofx/gistic_r:v1"
        memory: "8 GB"
    }
}

task run_gistic2 {
    #Inputs defined here
    File seg_file
    File markers_file
    File refgene_file
    String prefix

    Float? cap=1.5
    Float? broad_length_cutoff=0.5
    Int? remove_X=1
    Float? conf_level=0.90
    Int? join_segment_size=4
    Int? arm_peel=1
    Int? max_sample_segs=2500
    Int? do_gene_gistic=1
    String? gene_collapse_method="extreme"
    Int? memoryGB=8
    Int? preemptible=0

    parameter_meta {

        seg_file: "A six-column tab-delimited file containing the segmented data for all tumor/normal pairs in the pair set."
        markers_file: "A three-column tab-delimited file identifying the names and positions all markers."
        refgene_file: "Contains information about the location of genes and cytobands on a given build of the genome.  These files are created in MATLAB."
        cap: "Minimum and maximum cap values on analyzed data. Regions with a log2 ratio greater than the cap are set to the cap value; regions with a log2 ratio less than -cap value are set to -cap. (DEFAULT=1.5)"
        broad_length_cutoff: "Threshold used to distinguish broad from focal events, given in units of fraction of chromsome arm.  (Recommended: 0.7)"
        remove_X: "0/1 flag indicating whether to remove data from the X chromosome before analysis.  (Recommended: 0)"
        conf_level: "Confidence level used to calculate region containing the driver.  (Recommended: 0.99)"
        join_segment_size: "Smallest number of markers to allow in segements from the segemented data.  Segements that contain a number of markers less than or equal to this number are joined to the adjacent segement, closest in copy number. (Recommended: 4)"
        arm_peel: "0/1 flag indicating whether to perform arm-level peel off, wich helps spearte peaks and clean up noise.  (Recommended: 1)"
        max_sample_segs: "Maximum number of segements allowed for a smaple in the input data.    Samples with more segments than this are excluded from the analysis. (Recommended: 2000)"
        do_gene_gistic: "0/1 flag indicating tht the gene GISTIC aglrithm should be used to calculate signficance of deletions at the gene level instead of a marker level. (Recommended: 1)"
        gene_collapse_method: "Method for reducing marker-level copy number data to the gene-level copy number data in the gene tables. Markers contained in the gene are used when available, otherwise the flanking marker or markers are used. Allowed values are mean, median, min, max or extreme. The extreme method chooses whichever of min or max is furthest from diploid. (Recommended: extreme)"
        memoryGB: "Integer value specifying the minimum memory requirements (in GB) for the virtual machine running the GISTIC2 task."
        preemptible: "Integer value specifying the maximum number of times Cromwell should request a preemptible machine for this task before defaulting back to a non-preemptible one."
    }

    command {

    MCR_ROOT=/usr/local/MATLAB/MATLAB_Compiler_Runtime
    MCR_VER=v83

    echo Setting Matlab MCR root to $MCR_ROOT

    ## set up environment variables
    LD_LIBRARY_PATH=$MCR_ROOT/$MCR_VER/runtime/glnxa64:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=$MCR_ROOT/$MCR_VER/bin/glnxa64:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=$MCR_ROOT/$MCR_VER/sys/os/glnxa64:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH
    XAPPLRESDIR=$MCR_ROOT/$MCR_VER/MATLAB_Component_Runtime/v83/X11/app-defaults
    export XAPPLRESDIR

    mkdir gistic_out
    gp_gistic2_from_seg -b gistic_out -seg ${seg_file} -mk ${markers_file} -refgene ${refgene_file} -cap ${cap} -genegistic 1 -smallmem 1 -broad 1 -brlen ${broad_length_cutoff} -rx ${remove_X} -conf ${conf_level} -js ${join_segment_size} -armpeel ${arm_peel} -maxseg ${max_sample_segs} -savegene ${do_gene_gistic} -gcm ${gene_collapse_method}
    tar zcvf ${prefix}_gistic_out.tar.gz gistic_out
    }

    output {
        File gene_level_cnv="gistic_out/all_thresholded.by_genes.txt"
        File gistic_out_tar="${prefix}_gistic_out.tar.gz"
    }

    runtime {
        docker: "genepattern/docker-gistic:0.12"
        slurm_docker: "genepattern/docker-gistic:0.12"
        memory: "${memoryGB} GB"
        preemptible: "${preemptible}"
    }

    meta {
         author : "Xihao Hu"
         email : "huxihao@gmail.com"
    }

}

workflow Gistic {
    Array[File] seq_files
    String prefix

    call seq_to_gistic{
        input:
           seq_files=seq_files,
    }
    call run_gistic2{
        input:

            prefix=prefix,
            markers_file=seq_to_gistic.markers_file,
            seg_file=seq_to_gistic.seq_file
    }
}
