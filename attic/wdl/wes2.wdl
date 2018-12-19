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

    # Pipeline for HLA-typing
    call polysolver_wdl.strip_chr_substring {
       input: input_bam = normal_align_recalibrated_bam,
              memory=4,
              num_cpu=8,
              num_preempt=0,
              disk_space=60,
              boot_disk_gb=20
    }
    call polysolver_wdl.hla_type {
       input: bam=strip_chr_substring.output_bam,
              bam_index=strip_chr_substring.output_bam_index,
              build="hg38",
              memory=4,
              num_cpu=8,
              num_preempt=0,
              disk_space=30,
              boot_disk_gb=10,
              prefix=prefix
    }

    # Pipeline for neo-antigen prediction
    call oncotator_wdl.liftover_vcf {
       input: input_vcf = run_mutect2_unfiltered_vcf,
              liftover_chain = liftover_hg38_hg19_chain,
              reference_fasta = liftover_hg19_fasta,
              reference_dict = liftover_hg19_dict,
              prefix = prefix
    }
    call oncotator_wdl.oncotate_m2 {
       input: m2_vcf = liftover_vcf.valid_vcf,
              case_id = prefix+'.tumor',
              control_id = prefix+'.normal',
              onco_ds_tar_gz = oncotate_ds_tar_gz,
              default_config_file = oncotate_config_file
    }
    call neoantigen_wdl.Neoantigen as neoantigen {
       input: hg19DBTarBall = neoantigen_hg19DBTarBall,
              hlaFile = hla_type.winners_file,
              mutFile = oncotate_m2.oncotated_m2_maf,
              indelFile = oncotate_m2.oncotated_m2_maf,
              id = prefix,
              preemptible = 0,
              rankVersion = "2",
              binderThreshold = "all",
              hlaCStatus = "1" 
    }

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
        # HLA-typing
        File hla_type_winners_file = hla_type.winners_file
        Array[Array[String]] hla_type_winners = hla_type.winners

        # Neoantigen
        File oncotate_m2_maf = oncotate_m2.oncotated_m2_maf
        File neoantigen_combinedFullpeptideFile = neoantigen.combinedFullpeptideFile
        File neoantigen_combinedAllBindersCleanFile = neoantigen.combinedAllBindersCleanFile
        File neoantigen_outFilterInconsistentFile = neoantigen.outFilterInconsistentFile

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
