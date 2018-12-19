import "../cidc-pipelines/wdl/wes/wes2_partial.wdl" as wes2_wdl
import "../cidc-pipelines/wdl/wes/gistic.wdl" as gistic

workflow run_multi_wes2 {
    # input bams must be name sorted
    String cohort
    Array[Map[String,File]]  patients
    scatter(patient in patients) {
        call wes2_wdl.run_wes2 as run_wes2 {
            input:
               normal_align_recalibrated_bam = patient['normal_bam'],
               normal_align_recalibrated_bam_index = patient['normal_bai'],
               run_mutect2_unfiltered_vcf= patient['vcf'],
               tumor_align_recalibrated_bam = patient['tumor_bam'],
               tumor_align_recalibrated_bam_index = patient['normal_bam'],
               prefix = patient['prefix']
        }
    }
    call gistic.Gistic as run_gistic{
        input:
            seq_files = run_wes2.sequenza_gistic_seg,
            prefix = cohort
    }
}

