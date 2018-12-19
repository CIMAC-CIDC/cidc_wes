import "../cidc-pipelines/wdl/utilities/samtosortedbam.wdl" as stb
import "../cidc-pipelines/wdl/utilities/idxstats.wdl" as idx

task bwamem {
    File fastq1
    File? fastq2
    Array[File] index
    String? prefix = "output"

    Int? gzipped = 1 # 1 for yes or 0 for no expect gzipped fastq's by default
    Int? output_all_secondary = 0 # if 1 for yes will add the -a flag
    Int? sort_threads = 8
    Int? sort_memory_mb = 768
    Int? complevel = 9

    # Required options
    Int memory
    Int disk_space
    Int num_cpu
    Int? num_preempt = 0
    Int? boot_disk_gb = 20

    String? docker = "vacation/bwasam:0.7.15"
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
    parameter_meta {
       fastq1: "fastq file to analyze"
    }
    command {
       mkdir mytemp
       INDEX=$(echo $(ls ${index[1]} | head -n 1) | rev |  cut -d '.' -f 2- | rev)
       ALLOPTION=""
       if [ ${output_all_secondary}==1 ]
       then
           ALLOPTION="-a"
       fi
       if [ -f ${fastq2} ]
       then
           if [ ${gzipped}==1 ]
           then
               bwa mem -t ${num_cpu} \
                   $ALLOPTION \
                   -M \
                   $INDEX \
                   <(gunzip -c ${fastq1}) \
                   <(gunzip -c ${fastq2}) \
                   | samtools view -Sb -u - \
                   | samtools sort -m ${sort_memory_mb}M -T $(pwd)/mytemp/mysort -@ ${sort_threads} -l ${complevel} \
                                   - -o "${prefix}.sorted.bam"
           else
               bwa mem -t ${num_cpu} \
                   $ALLOPTION \
                   -M \
                   $INDEX \
                   ${fastq1} \
                   ${fastq2} \
                   | samtools view -Sb -u - \
                   | samtools sort -m ${sort_memory_mb}M -T $(pwd)/mytemp/mysort -@ ${sort_threads} -l ${complevel} \
                                   - -o "${prefix}.sorted.bam"
           fi
       else
           if [ ${gzipped}==1 ]
           then
               bwa mem -t ${num_cpu} \
                   $ALLOPTION \
                   -M \
                   $INDEX \
                   <(gunzip -c ${fastq1}) \
                   | samtools view -Sb -u - \
                   | samtools sort -m ${sort_memory_mb}M -T $(pwd)/mytemp/mysort -@ ${sort_threads} -l ${complevel} \
                                   - -o "${prefix}.sorted.bam"
           else
               bwa mem -t ${num_cpu} \
                   $ALLOPTION \
                   -M \
                   $INDEX \
                   ${fastq1} \
                   | samtools view -Sb -u - \
                   | samtools sort -m ${sort_memory_mb}M -T $(pwd)/mytemp/mysort -@ ${sort_threads} -l ${complevel} \
                                   - -o "${prefix}.sorted.bam"
           fi
       fi
       samtools index "${prefix}.sorted.bam"
       rm -r mytemp
    }
    output {
       File bam = "${prefix}.sorted.bam"
       File bam_index = "${prefix}.sorted.bam.bai"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN]"
    }
}

workflow run_bwamem {
    String? prefix = "output"

    Int? gzipped = 1 # 1 for yes or 0 for no expect gzipped fastq's by default
    Int? output_all_secondary = 0 # if 1 for yes will add the -a flag
    Int? complevel = 9

    call bwamem {
        input:
            prefix = prefix,
            complevel = complevel,
            output_all_secondary = output_all_secondary,
            gzipped = gzipped
    }
    call idx.idxstats as idxstats {
       input: bam=bwamem.bam,
              bam_index=bwamem.bam_index
    }
    output {
        File bam = bwamem.bam
        File bam_index = bwamem.bam_index
        Array[Array[String]] idxstats_table = idxstats.idxstats_table
        Float fraction_aligned_reads = idxstats.fraction_aligned_reads
        Int total_reads = idxstats.total_reads
        Int aligned_reads = idxstats.aligned_reads
        Int unaligned_reads = idxstats.unaligned_reads
    }
}

