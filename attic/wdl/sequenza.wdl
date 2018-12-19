task sequenza {
  String sampleName
  File normalBam
  File tumorBam
  File ReferenceFasta
  File ReferenceGcWig
  String? assembly = "hg38"
  String? docker = "gcr.io/cidc-biofx/sequenza:v1"
  String? num_cpu = "8"
  String? memory = "16"

  command <<<

  set -x
  sequenza-pipeline \
      --sample-id ${sampleName} \
      --normal-bam ${normalBam} \
      --tumor-bam ${tumorBam} \
      --reference-gz ${ReferenceFasta} \
      --gc_wig ${ReferenceGcWig}

  cp seqz/${sampleName}_bin50.seqz.gz ${sampleName}.bin50.txt.gz

  cat <<EOT >> run_seq.R

library(sequenza)

test <- sequenza.extract("${sampleName}.bin50.txt.gz", assembly = "${assembly}")
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = "sample", out.dir="OUT")

mut.tab <- read.table("OUT/sample_mutations.txt", header=TRUE)
seg.res <- read.table("OUT/sample_segments.txt", header=TRUE)

## fix chr list. Remark: do not use '$'

chrs = names(test[['mutations']])
mut.tab[['chromosome']] = factor(as.character(mut.tab[['chromosome']]), levels=chrs)
seg.res[['chromosome']] = factor(as.character(seg.res[['chromosome']]), levels=chrs)

tsv <- sequenza:::sequenza2PyClone(mut.tab, seg.res, "${sampleName}", norm.cn = 2)
tsv <- tsv[tsv[,'major_cn'] != 0 ,] # Required for Pyclone
write.table(tsv, file = "${sampleName}.pyclone.tsv", sep='\t', row.names=FALSE, quote=FALSE)

EOT

  Rscript run_seq.R
  tar zcvf ${sampleName}.seqz_fit.tar.gz OUT/
  cp OUT/sample_segments.txt ${sampleName}_segments.txt

  >>>

  runtime {
    docker: "${docker}"
    slurm_docker: "${docker}"
    memory: "${memory} GB"
    slurm_memory: "${memory}G"
    cpu: "${num_cpu}"
    disks: "local-disk 100 HDD"
    preemptible: "0"
    zones: "us-east1-b us-east1-c us-east1-d"
  }

  output {
    File bin_seqz_gz = "${sampleName}.bin50.txt.gz"
    File seqz_fit_tar = "${sampleName}.seqz_fit.tar.gz"
    File pyclone_tsv = "${sampleName}.pyclone.tsv"
    File gistic_seg = "${sampleName}_segments.txt"
  }
}

workflow Sequenza {
  String sampleName
  File normalBam
  File tumorBam
  File ReferenceFasta
  File ReferenceGcWig

  call sequenza {
    input: sampleName=sampleName,
           normalBam=normalBam, tumorBam=tumorBam,
           ReferenceFasta=ReferenceFasta,
           ReferenceGcWig=ReferenceGcWig
  }
}
