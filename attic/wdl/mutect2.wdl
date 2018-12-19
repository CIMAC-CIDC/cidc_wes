import "../cidc-pipelines/wdl/utilities/split_bam_sequences.wdl" as split_bam_sequences_wdl
import "../cidc-pipelines/wdl/utilities/merge_vcfs.wdl" as merge_vcfs_wdl
import "../cidc-pipelines/wdl/utilities/sort_compress_vcf.wdl" as sort_compress_vcf_wdl

# adapted from https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/blob/master/mutect2.wdl
#
# License: BSD 3-Clause License

#  Run Mutect 2 on a single tumor-normal pair or on a single tumor sample.
#
#  Description of inputs
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  tumor_bam, tumor_bam_index: self-explanatory
#  normal_bam, normal_bam_index: self-explanatory
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_index: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_index: vcf of common variants with allele frequencies fo calculating contamination
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  is_run_oncotator: if true, annotate the mutect2 VCFs using oncotator (to produce a TCGA MAF).  Important:  This requires a docker image and should
#   not be run in environments where docker is unavailable (e.g. SGE cluster on a Broad on-prem VM).  Access to docker
#   hub is also required, since the task will download a public docker image.
#
#
workflow run_mutect2 {
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  File normal_bam
  File normal_bam_index
  String prefix
  Int num_threads
  File? pon
  File? pon_index
  File gnomad
  File gnomad_index
  File? variants_for_contamination
  File? variants_for_contamination_index
  Boolean is_run_orientation_bias_filter

  Array[String]? artifact_modes
  String? m2_extra_args
  String? m2_extra_filtering_args
  String? sequencing_center
  String? sequence_source
  File? default_config_file
  Boolean is_bamOut = false
  Boolean read_groups_set = false

  Int? preemptible_attempts
  String gatk_docker = "vacation/gatk:4.0.3.0"

  call split_bam_sequences_wdl.get_sequence_list {
      input:
          input_bam = tumor_bam
  }
  scatter (sequence_string in get_sequence_list.sequence_list) {
      call split_bam_sequences_wdl.get_bam_sequence as get_tumor {
          input:
              input_bam = tumor_bam,
              input_bam_index = tumor_bam_index,
              sequence_string = sequence_string
      }
      call split_bam_sequences_wdl.get_bam_sequence as get_normal {
          input:
              input_bam = normal_bam,
              input_bam_index = normal_bam_index,
              sequence_string = sequence_string
      }
  }
  ## Run the operations for everything except the unaligned reads
  scatter (i in range(length(get_tumor.output_bam_and_index))) {
    call mutect2 {
      input: 
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        tumor_bam = get_tumor.output_bam_and_index[i]["bam"],
        tumor_bam_index = get_tumor.output_bam_and_index[i]["index"],
        normal_bam = get_normal.output_bam_and_index[i]["bam"],
        normal_bam_index = get_normal.output_bam_and_index[i]["index"],
        intervals = get_tumor.output_bam_and_index[i]["interval"],
        pon = pon,
        pon_index = pon_index,
        gnomad = gnomad,
        gnomad_index = gnomad_index,
        preemptible_attempts = preemptible_attempts,
        #dbsnp = dbsnp,
        #dbsnp_index = dbsnp_index,
        m2_extra_args = m2_extra_args,
        is_bamOut = is_bamOut,
        preemptible_attempts = preemptible_attempts,
        gatk_docker = gatk_docker,
        num_threads = num_threads,
        read_groups_set = read_groups_set
    }
  }
  call merge_vcfs_wdl.merge_vcfs as merge_vcfs {
      input: 
          input_vcfs = mutect2.output_vcf,
          prefix = prefix
  }
  if (is_run_orientation_bias_filter) {
      call CollectSequencingArtifactMetrics {
        input:
            preemptible_attempts = preemptible_attempts,
            gatk_docker = gatk_docker,
            tumor_bam = tumor_bam,
            tumor_bam_index = tumor_bam_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index
      }
  }
  call Filter {
    input:
      unfiltered_vcf = merge_vcfs.output_vcf,
      unfiltered_vcf_index = merge_vcfs.output_vcf_index,
      gatk_docker = gatk_docker,
      preemptible_attempts = preemptible_attempts,
      pre_adapter_metrics = CollectSequencingArtifactMetrics.pre_adapter_metrics,
      tumor_bam = tumor_bam,
      tumor_bam_index = tumor_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      artifact_modes = artifact_modes,
      variants_for_contamination = variants_for_contamination,
      variants_for_contamination_index = variants_for_contamination_index,
      m2_extra_filtering_args = m2_extra_filtering_args
  }
  call sort_compress_vcf_wdl.sort_compress_vcf as sort_compress_vcf {
      input:
          input_vcf = Filter.filtered_vcf,
          prefix = prefix+'.filtered'
  }
  output {
        File unfiltered_vcf = merge_vcfs.output_vcf
        File unfiltered_vcf_index = merge_vcfs.output_vcf_index
        File filtered_vcf = sort_compress_vcf.output_vcf
        File filtered_vcf_index = sort_compress_vcf.output_vcf_index
        File contamination_table = Filter.contamination_table

        # select_first() fails if nothing resolves to non-null, so putting in "null" for now.
        File? preadapter_detail_metrics = select_first([CollectSequencingArtifactMetrics.pre_adapter_metrics, "null"])
        #Array[File]? bamout = mutect2.output_bamOut
  }
}

task mutect2 {
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  File? normal_bam
  File? normal_bam_index
  File? pon
  File? pon_index
  File? gnomad
  File? gnomad_index
  String? m2_extra_args
  Boolean? is_bamOut
  Int num_threads
  Boolean? read_groups_set

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command <<<

   # We need to create these files regardless, even if they stay empty
   touch bamout.bam
   echo "" > normal_name.txt

  if [[ "${read_groups_set}" != "true" ]]
  then
    gatk --java-options "-Xmx4g" AddOrReplaceReadGroups -I ${tumor_bam} \
      -O tumor_retag.bam \
      --RGID tumor \
      --RGLB tumor \
      --RGPL tumor \
      --RGPU tumor \
      --RGSM tumor
    samtools index tumor_retag.bam

    gatk --java-options "-Xmx4g" AddOrReplaceReadGroups -I ${normal_bam} \
      -O normal_retag.bam \
      --RGID normal \
      --RGLB normal \
      --RGPL normal \
      --RGPU normal \
      --RGSM normal
    samtools index normal_retag.bam
  gatk --java-options "-Xmx4g" Mutect2 \
    -R ${ref_fasta} \
    --native-pair-hmm-threads ${num_threads} \
    -I tumor_retag.bam \
    --tumor-sample tumor \
    -I normal_retag.bam \
    --normal-sample normal \
    ${"--germline-resource " + gnomad} \
    ${"-pon " + gnomad} \
    ${"-L " + intervals} \
    -O "output.vcf" \
    ${true='--bam-output bamout.bam' false='' is_bamOut} \
    ${m2_extra_args}
  else
  gatk --java-options "-Xmx4g" Mutect2 \
    -R ${ref_fasta} \
    --native-pair-hmm-threads ${num_threads} \
    -I ${tumor_bam} \
    --tumor-sample tumor \
    -I ${normal_bam} \
    --normal-sample normal \
    ${"--germline-resource " + gnomad} \
    ${"-pon " + gnomad} \
    ${"-L " + intervals} \
    -O "output.vcf" \
    ${true='--bam-output bamout.bam' false='' is_bamOut} \
    ${m2_extra_args}
  fi
  
  >>>

  runtime {
        docker: "${gatk_docker}"
        slurm_docker: "${gatk_docker}"
        memory: select_first([mem, 12]) + " GB"
        slurm_memory: select_first([mem, 12]) + "G"
        disks: "local-disk " + select_first([disk_space_gb, 200]) + " HDD"
        bootDiskSizeGb: "50"
        preemptible: select_first([preemptible_attempts, 2])
        zones: "us-east1-b us-east1-c us-east1-d"
        cpu: "${num_threads}"
  }

  output {
    File output_vcf = "output.vcf"
    File output_bamOut = "bamout.bam"
  }
  meta {
        maintainer: "Jason L Weirather"
        originalsource: "https://github.com/gatk-workflows/gatk4-somatic-snvs-indels/blob/master/mutect2.wdl"
        license: "BSD 3-Clause License"
        citation: "Cibulskis K, Lawrence MS, Carter SL, Sivachenko A, Jaffe D, Sougnez C, Gabriel S, Meyerson M, Lander ES, Getz G. Sensitive detection of somatic point mutations  in impure and heterogeneous cancer samples. Nat Biotechnol. 2013 Mar;31(3):213-9. doi: 10.1038/nbt.2514. Epub 2013 Feb 10. PubMed PMID: 23396013; PubMed Central PMCID: PMC3833702."
  }
}

task CollectSequencingArtifactMetrics {
  File tumor_bam
  File tumor_bam_index
  File ref_fasta
  File ref_fasta_index

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command {
        set -e
        gatk --java-options "-Xmx4g" CollectSequencingArtifactMetrics -I ${tumor_bam} -O "gatk" -R ${ref_fasta} --VALIDATION_STRINGENCY LENIENT
  }

  runtime {
        docker: "${gatk_docker}"
        slurm_docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        slurm_memory: select_first([mem, 5]) + "G"
        disks: "local-disk " + select_first([disk_space_gb, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
        zones: "us-east1-b us-east1-c us-east1-d"
  }

  output {
    File pre_adapter_metrics = "gatk.pre_adapter_detail_metrics"
  }
}

task Filter {
  File unfiltered_vcf
  File unfiltered_vcf_index
  String filtered_vcf_name = basename(unfiltered_vcf, ".sorted.vcf.gz") + "-filtered.vcf"
  File? pre_adapter_metrics
  File? tumor_bam
  File? tumor_bam_index
  File? ref_fasta
  File? ref_fasta_index
  Array[String]? artifact_modes
  File? variants_for_contamination
  File? variants_for_contamination_index
  String? m2_extra_filtering_args

  # Runtime parameters
  Int? mem
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb

  command {
    set -e

    # Use GATK Jar override if specified

    touch contamination.table
    if [[ "${variants_for_contamination}" == *.vcf ]]; then
        gatk --java-options "-Xmx4g" GetPileupSummaries -I ${tumor_bam} -V ${variants_for_contamination} -O pileups.table
        gatk --java-options "-Xmx4g" CalculateContamination -I pileups.table -O contamination.table
        contamination_cmd="--contamination-table contamination.table"
    fi

    gatk --java-options "-Xmx4g" FilterMutectCalls -V ${unfiltered_vcf} \
         -O filtered.vcf $contamination_cmd \
         ${m2_extra_filtering_args}

    # FilterByOrientationBias must come after all of the other filtering.
    if [[ ! -z "${pre_adapter_metrics}" ]]; then
        gatk --java-options "-Xmx4g" FilterByOrientationBias -AM ${sep=" -AM " artifact_modes} \
            -V filtered.vcf -P ${pre_adapter_metrics} --output ${filtered_vcf_name}
    else
        mv filtered.vcf ${filtered_vcf_name}
        mv filtered.vcf.idx "${filtered_vcf_name}.idx"
    fi
  }

  runtime {
        docker: "${gatk_docker}"
        slurm_docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        slurm_memory: select_first([mem, 5]) + "G"
        disks: "local-disk " + select_first([disk_space_gb, 200]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
        zones: "us-east1-b us-east1-c us-east1-d"
  }

  output {
    File filtered_vcf = "${filtered_vcf_name}"
    File filtered_vcf_index = "${filtered_vcf_name}.idx"
    File contamination_table = "contamination.table"
  }
}
