import "../cidc-pipelines/wdl/wes/oncotator.wdl" as oncotator_wdl
import "../cidc-pipelines/wdl/wes/neoantigen_mt.wdl" as neoantigen_wdl
import "../cidc-pipelines/wdl/wes/polysolver.wdl" as polysolver_wdl
import "../cidc-pipelines/wdl/wes/facets.wdl" as facets_wdl
import "../cidc-pipelines/wdl/wes/sequenza.wdl" as sequenza_wdl
import "../cidc-pipelines/wdl/wes/pyclone.wdl" as pyclone_wdl


workflow run_wes2 {
    String prefix

    # WES saved files
    File tumor_align_recalibrated_bam
    File tumor_align_recalibrated_bam_index
    File normal_align_recalibrated_bam
    File normal_align_recalibrated_bam_index
    File run_mutect2_unfiltered_vcf
    # Additional files
    File facets_vcftar
    File liftover_hg38_hg19_chain
    File liftover_hg19_fasta
    File liftover_hg19_dict
    File oncotate_ds_tar_gz
    File oncotate_config_file
    File neoantigen_hg19DBTarBall
    File reference_fasta
    File reference_gc_wiggle

    # Pipeline for tumor purity
    call facets_wdl.facets {
       input: tumorBam = tumor_align_recalibrated_bam,
              tumorBamIndex = tumor_align_recalibrated_bam_index,
              normalBam = normal_align_recalibrated_bam,
              normalBamIndex = normal_align_recalibrated_bam_index,
              vcftar = facets_vcftar,
              pairName = prefix
    }

    # Pipeline for clonality
    call sequenza_wdl.sequenza {
        input: sampleName=prefix,
               normalBam=normal_align_recalibrated_bam,
               tumorBam=tumor_align_recalibrated_bam,
               ReferenceFasta=reference_fasta,
               ReferenceGcWig=reference_gc_wiggle
    }
    call pyclone_wdl.pyclone {
       input: prefix=prefix,
              pyclone_tsv=sequenza.pyclone_tsv,
              memory=8
    }
    output {
        # Facets
        File facets_facetsOutput = facets.facetsOutput
        Array[Array[String]] facets_facetsEstimation = facets.facetsEstimation
        File facets_facetsIterations = facets.facetsIterations
        File facets_pileup = facets.pileup
        File facets_cncf = facets.cncf
        File facets_plot = facets.plot

        # Sequenza
        File sequenza_bin_seqz_gz = sequenza.bin_seqz_gz
        File sequenza_seqz_fit_tar = sequenza.seqz_fit_tar
        File sequenza_pyclone_tsv = sequenza.pyclone_tsv
        File sequenza_gistic_seg = sequenza.gistic_seg

        # Pyclone
        File pyclone_pyclone_tar = pyclone.pyclone_tar

    }
    meta {
        maintainer: "Xihao Sherlock Hu"
    }
}
