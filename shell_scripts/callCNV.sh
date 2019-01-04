#!/bin/sh
cohort=$1
patient=$2
TARGET_BED="/cluster/jxfu/proj/CIDC/Sentieon/data/Reference/$3"
TUMOR_BAM="/cluster/jxfu/proj/CIDC/Sentieon/data/mutation/${cohort}/$patient/tumor_recalibrated.bam"
NUMBER_THREADS=8
REFERENCE="/cluster/asahu/mutation_calling/MDAnderson/ref/Homo_sapiens_assembly38.fasta"
PADDING=0
PON="/cluster/jxfu/proj/CIDC/Sentieon/data/PoN/${cohort}"
OUT_CNV_folder="/cluster/jxfu/proj/CIDC/Sentieon/data/CNV/${cohort}"
OUT_CNV="/cluster/jxfu/proj/CIDC/Sentieon/data/CNV/${cohort}/result/$patient"
OUT_Gistic="/cluster/jxfu/proj/CIDC/Sentieon/data/CNV/${cohort}/comp_sequanza/${patient}_segments.txt"
release_dir="/cluster/jxfu/proj/CIDC/Sentieon/release/sentieon-genomics-201808.01"
logfile="pipe_log/cnv_${cohort}_${patient}.log"
#0.Setup
export SENTIEON_LICENSE=172.24.216.24:8990
mkdir -p $OUT_CNV_folder/result
mkdir -p $OUT_CNV_folder/comp_sequanza

exec >$logfile 2>&1

$release_dir/bin/sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
-i $TUMOR_BAM --algo CNV \
--target $TARGET_BED \
--target_padding $PADDING --pon $PON \
$OUT_CNV

## Convert output to gistic input
(echo -e "\"chromosome\"\t\"start.pos\"\t\"end.pos\"\t\"N.BAF\"\t\"depth.ratio\"" ; awk '($2 !~ "_") && ( FNR > 1) {print $2"\t"$3"\t"$4"\t"$5"\t"$6}' $OUT_CNV ) &> $OUT_Gistic
