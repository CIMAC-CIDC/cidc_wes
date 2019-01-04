#!/bin/sh
# *******************************************
# Script to perform TN seq variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************
cohort=$1
patient=$2
tumor_suffix=$3
cohort=${cohort%$tumor_suffix}
#fastq_folder="/cluster/cidc/Immunotherapy/WES/${cohort}/merge_fastq/" For HNSC only
fastq_folder="/cluster/cidc/Immunotherapy/WES/${cohort}/fastq/"
tumor_fastq_1="${patient}.T${tumor_suffix}.1.fq.gz"
tumor_fastq_2="${patient}.T${tumor_suffix}.2.fq.gz"
tumor_sample="${patient}.T${tumor_suffix}"
tumor_group="${patient}.T${tumor_suffix}"
normal_fastq_1="${patient}.N.1.fq.gz"
normal_fastq_2="${patient}.N.2.fq.gz " #If using Illumina paired data
normal_sample="${patient}.N"
normal_group="${patient}.N"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/cluster/asahu/mutation_calling/MDAnderson/ref/Homo_sapiens_assembly38.fasta"
dbsnp="/cluster/asahu/mutation_calling/MDAnderson/ref/Homo_sapiens_assembly38.dbsnp138.vcf"
known_Mills_indels="/cluster/asahu/mutation_calling/MDAnderson/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_1000G_indels="/cluster/asahu/mutation_calling/script/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# Update with the location of the Sentieon software package and license file
release_dir="/cluster/jxfu/proj/CIDC/Sentieon/release/sentieon-genomics-201808.01"
export SENTIEON_LICENSE=172.24.216.24:8990

# Other settings
nt=8 #number of threads to use in computation
workdir="/cluster/jxfu/proj/CIDC/Sentieon/data/mutation/${cohort}/${patient}${tumor_suffix}" #Determine where the output files will besored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000
( $release_dir/bin/bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$tumor_fastq_1 $fastq_folder/$tumor_fastq_2 || echo -n 'error' ) | $release_dir/bin/sentieon util sort -o tumor_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $release_dir/bin/bwa mem -M -R "@RG\tID:$normal_group\tSM:$normal_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$normal_fastq_1 $fastq_folder/$normal_fastq_2 || echo -n 'error' ) | $release_dir/bin/sentieon util sort -o normal_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_sorted.bam --algo MeanQualityByCycle tumor_mq_metrics.txt --algo QualDistribution tumor_qd_metrics.txt --algo GCBias --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
#$release_dir/bin/sentieon plot metrics -o tumor_metrics-report.pdf gc=tumor_gc_metrics.txt qd=tumor_qd_metrics.txt mq=tumor_mq_metrics.txt isize=tumor_is_metrics.txt
# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_sorted.bam --algo MeanQualityByCycle normal_mq_metrics.txt --algo QualDistribution normal_qd_metrics.txt --algo GCBias --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt
#$release_dir/bin/sentieon plot metrics -o normal_metrics-report.pdf gc=normal_gc_metrics.txt qd=normal_qd_metrics.txt mq=normal_mq_metrics.txt isize=normal_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo LocusCollector --fun score_info tumor_score.txt
$release_dir/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo Dedup --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam 
# ******************************************
# 3b. Remove Duplicate Reads for normal sample
# ******************************************
$release_dir/bin/sentieon driver -t $nt -i normal_sorted.bam --algo LocusCollector --fun score_info normal_score.txt
$release_dir/bin/sentieon driver -t $nt -i normal_sorted.bam --algo Dedup --rmdup --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam 

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tumor_realigned.bam
# ******************************************
# 4b. Indel realigner for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels normal_realigned.bam

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table --algo ReadWriter tumor_recalibrated.bam
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$release_dir/bin/sentieon plot bqsr -o tumor_recal_plots.pdf tumor_recal.csv
# ******************************************
# 5b. Base recalibration for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table --algo ReadWriter normal_recalibrated.bam

$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam -q normal_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table.post
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
$release_dir/bin/sentieon plot bqsr -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 6. Corealignment of tumor and normal
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tn_corealigned.bam

# ******************************************
# 7. Somatic Variant Calling
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tn_corealigned.bam --algo TNsnv --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp --call_stats_out output-call.stats output-tnsnv.vcf.gz
vcftools --gzvcf output-tnsnv.vcf.gz  --remove-filtered-all --recode --stdout > output-tnsnv.filter.vcf

$release_dir/bin/sentieon driver -r $fasta -t $nt -i tn_corealigned.bam --algo TNhaplotyper --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp output-tnhaplotyper.vcf.gz
vcftools --gzvcf output-tnhaplotyper.vcf.gz  --remove-filtered-all --recode --stdout > output-tnhaplotyper.filter.vcf

$release_dir/bin/sentieon driver -r $fasta -t $nt -i tn_corealigned.bam --algo TNscope --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp output-tnscope.vcf.gz
vcftools --gzvcf output-tnscope.vcf.gz  --remove-filtered-all --recode --stdout > output-tnscope.filter.vcf

