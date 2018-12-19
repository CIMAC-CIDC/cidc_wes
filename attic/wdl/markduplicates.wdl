task markduplicates {
    File input_bam
    String? prefix = "markdup"
    Int? OPTICAL_DUPLICATE_PIXEL_DISTANCE = 100
    String? ASSUME_SORT_ORDER = "coordinate"

    # Required options
    Int? memory = 32
    Int? disk_space = 100
    Int? num_cpu = 2
    Int? num_preempt = 0
    Int? boot_disk_gb = 50

    String? docker = "broadinstitute/gatk:4.0.1.2"

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
       input_bam: "bam in the proper format for processing by picard markduplicates suchas bwa mem -M"
    }
    command <<<
       mkdir mytemp
       gatk  --java-options "-Xmx16g" MarkDuplicates \
           --TMP_DIR $(pwd)/mytemp \
           -I ${input_bam} \
           -O ${prefix}.markduplicates.bam \
           --METRICS_FILE $(pwd)/${prefix}.markduplicates.bam.txt \
           --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${OPTICAL_DUPLICATE_PIXEL_DISTANCE} \
           --ASSUME_SORT_ORDER ${ASSUME_SORT_ORDER} \
           --VALIDATION_STRINGENCY LENIENT 2> /dev/null
        cat ${prefix}.markduplicates.bam.txt | grep -A 2 -e '^## METRICS' | tail -n 2 > myfile1
        awk '{ 
            for (i=1; i<=NF; i++)  {
                a[NR,i] = $i
            }
        }
       NF>p { p = NF }
       END {    
           for(j=1; j<=p; j++) {
               str=a[1,j]
               for(i=2; i<=NR; i++){
                   str=str"\t"a[i,j];
               }
           print str
       }
       }' myfile1 > myout
       rm -r mytemp
    >>>
    output {
       File output_bam = "${prefix}.markduplicates.bam"
       File metrics = "${prefix}.markduplicates.bam.txt"
       Map[String, String] summary = read_map("myout")       
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508."
    }
}

workflow run_markduplicates {
    call markduplicates {
    }
}
