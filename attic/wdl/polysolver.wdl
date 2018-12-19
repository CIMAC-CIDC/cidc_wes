# Input: Normal and Tumor WES
# Output: Polysolver HLA calls and/or HLA mutation calls
#
# The workflow defined in this docker assumes that your chromosome
# names are in the chr1, chr2 format rather than 1, 2 format which
# polysolver expects.  That is why there are extra steps to strip
# these names from the dockerfile, and a python script in our utilities
# does this from a dockerfile called polysolverhelp we build for this
# purpose.  If you are starting from BAM's aligned to chromosomes in the
# format polysolver expects you can omit these name striping steps.

task hla_type_mutations {
    File normal_bam
    File normal_bam_index
    File tumor_bam
    File tumor_bam_index
    File winners
    String build
    String? prefix = "output"

    # Required options
    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    String? docker = "vacation/polysolver:v4"

    parameter_meta {
       normal_bam: "Normal WES alignment"
       normal_bam_index: "Index of Normal WES alignment, presumably in the same location as the bam"
       tumor_bam: "Tumor WES alignment"
       tumor_bam_index: "Index of Tumor WES in same location as as the bam"
       winners: "File listing calls of HLA types"
       build: "Genome build hg38 or hg19"
       prefix: "Optional string prefix to give some of the output files."
    }

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
    command {
        echo $(pwd)
        export SAMTOOLS_DIR="/home/polysolver/binaries"
        bash /home/polysolver/scripts/shell_call_hla_mutations_from_type "${normal_bam}" "${tumor_bam}"  ${winners} ${build} STDFQ $(pwd) ${prefix}
    }
    output {
       File hla_mut_raw = "hla_mut.tar.gz"
       File hla_type_bams = "hla_type_bams.tar.gz"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "Shukla SA, Rooney MS, Rajasagi M, Tiao G, Dixon PM, Lawrence MS, Stevens J, Lane WJ, Dellagatta JL, Steelman S, Sougnez C, Cibulskis K, Kiezun A, Hacohen N,  Brusic V, Wu CJ, Getz G. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. Nat Biotechnol. 2015 Nov;33(11):1152-8. PubMed PMID: 26372948; PubMed Central PMCID: PMC4747795."
    }
}

task hla_type {
    File bam
    File bam_index
    String prefix
    String? race = "Unknown"
    Int? includeFreq = 1
    String build
    Int? insertCalc = 0

    # Required options
    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    String? docker = "vacation/polysolver:v4"

    parameter_meta {
       bam: "WES input alignment in bam format"
       bam_index: "index of bam in bai format and presumably the same location"
       race: "Unknown by default not sure how critical setting this is"
       includeFreq: "1 by default but need to ask Sachet about prefered defaults"
       build: "String hg38 or hg18 to describe the build of the genome"
       insertCalc: "0 by default, need to check with Sachet to make sure this is prefered"
    }

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
    command {
        echo $(pwd)
        export SAMTOOLS_DIR="/home/polysolver/binaries"
        bash /home/polysolver/scripts/shell_call_hla_type "${bam}" ${race} ${includeFreq} ${build} STDFQ ${insertCalc} $(pwd)
        mv winners.hla.txt ${prefix}.hla.txt
    }
    output {
       File winners_file = "${prefix}.hla.txt"
       Array[Array[String]] winners = read_tsv("${prefix}.hla.txt")
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "Shukla SA, Rooney MS, Rajasagi M, Tiao G, Dixon PM, Lawrence MS, Stevens J, Lane WJ, Dellagatta JL, Steelman S, Sougnez C, Cibulskis K, Kiezun A, Hacohen N,  Brusic V, Wu CJ, Getz G. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. Nat Biotechnol. 2015 Nov;33(11):1152-8. PubMed PMID: 26372948; PubMed Central PMCID: PMC4747795."
    }
}

task annotate_hla {
    File mutation_tar_gz
    String? prefix = "output"

    # Required options
    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    String? docker = "vacation/polysolver:v4"

    parameter_meta {
      mutation_tar_gz: "The tarball of mutation data that was produced by the hla mutation calling"
      prefix: "The same prefix used in mutation calling, by default we are using 'output'"
    }

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
    command {
        echo $(pwd)
        export SAMTOOLS_DIR="/home/polysolver/binaries"
        bash /home/polysolver/scripts/shell_annotate_hla_mutations ${prefix} ${mutation_tar_gz} $(pwd)
    }
    output {
       File mutect_ambiguous_annotated = "${prefix}.mutect.ambiguous.annotated"
       File mutect_filtered_annotated = "${prefix}.mutect.filtered.annotated"
       File mutect_filtered_nonsyn_annotated = "${prefix}.mutect.filtered.nonsyn.annotated"
       File mutect_filtered_syn_annotated = "${prefix}.mutect.filtered.syn.annotated"
       File mutect_unfiltered_annotated = "${prefix}.mutect.unfiltered.annotated"
       File strelka_indels_ambiguous_annotated = "${prefix}.strelka_indels.ambiguous.annotated"
       File strelka_indels_filtered_annotated = "${prefix}.strelka_indels.filtered.annotated"
       File strelka_indels_unfiltered_annotated = "${prefix}.strelka_indels.unfiltered.annotated"
    }
    meta {
        maintainer: "Jason L Weirather"
        citation: "Shukla SA, Rooney MS, Rajasagi M, Tiao G, Dixon PM, Lawrence MS, Stevens J, Lane WJ, Dellagatta JL, Steelman S, Sougnez C, Cibulskis K, Kiezun A, Hacohen N,  Brusic V, Wu CJ, Getz G. Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. Nat Biotechnol. 2015 Nov;33(11):1152-8. PubMed PMID: 26372948; PubMed Central PMCID: PMC4747795."
    }
}

task strip_chr_substring {
    File input_bam

    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    String? docker = "vacation/polysolverhelper:1.0.0"

    parameter_meta {
      input_bam: "The bam file which you will strip the 'chr' characters out of to make it compatible with polysolver"
    }

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
    command {
        echo $(pwd)
        python /remove_chr_substring_bam.py ${input_bam} output.bam
        samtools index output.bam
    }
    output {
       File output_bam = "output.bam"
       File output_bam_index = "output.bam.bai"
    }
}

workflow run_polysolver {
    File normal_bam
    File tumor_bam

    String build
    String? prefix = "output"

    Int memory
    Int disk_space
    Int num_cpu
    Int num_preempt
    Int boot_disk_gb

    call strip_chr_substring {
       input: input_bam = normal_bam,
              memory=memory,
              disk_space=disk_space,
              num_cpu=num_cpu,
              num_preempt=num_preempt,
              boot_disk_gb=boot_disk_gb
    }
    call strip_chr_substring as strip_chr_substring2 {
       input: input_bam = tumor_bam,
              memory=memory,
              disk_space=disk_space,
              num_cpu=num_cpu,
              num_preempt=num_preempt,
              boot_disk_gb=boot_disk_gb
    }
    call hla_type {
       input: bam=strip_chr_substring.output_bam,
              bam_index=strip_chr_substring.output_bam_index,
              build=build,
              memory=memory,
              disk_space=disk_space,
              num_cpu=num_cpu,
              num_preempt=num_preempt,
              boot_disk_gb=boot_disk_gb
    }
    call hla_type_mutations {
       input: normal_bam=strip_chr_substring.output_bam,
              normal_bam_index=strip_chr_substring.output_bam_index,
              tumor_bam=strip_chr_substring2.output_bam,
              tumor_bam_index=strip_chr_substring2.output_bam_index,
              winners=hla_type.winners_file,
              build=build,
              prefix=prefix,
              memory=memory,
              disk_space=disk_space,
              num_cpu=num_cpu,
              num_preempt=num_preempt,
              boot_disk_gb=boot_disk_gb
    }
    call annotate_hla {
       input: mutation_tar_gz = hla_type_mutations.hla_mut_raw,
              prefix = prefix,
              memory=memory,
              disk_space=disk_space,
              num_cpu=num_cpu,
              num_preempt=num_preempt,
              boot_disk_gb=boot_disk_gb
    }
}
