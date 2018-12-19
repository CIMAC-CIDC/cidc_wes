import "../cidc-pipelines/wdl/wes/wes.wdl" as wes_wdl

workflow run_multi_wes {
    # input bams must be name sorted
    Array[Map[String,File]]  fastqs
    scatter(fastq in fastqs) {
        call wes_wdl.run_wes as run_wes {
            input:
               tumor_fastq1 = fastq['tumor1'],
               tumor_fastq2 = fastq['tumor2'],
               normal_fastq1 = fastq['normal1'],
               normal_fastq2 = fastq['normal2'],
               prefix = fastq['prefix']
        }
    }
}

