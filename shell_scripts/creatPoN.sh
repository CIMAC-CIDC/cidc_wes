#!/bin/sh
cohort=$1
TARGET_BED="/cluster/jxfu/proj/CIDC/Sentieon/data/Reference/$2"
NUMBER_THREADS=8
REFERENCE="/cluster/asahu/mutation_calling/MDAnderson/ref/Homo_sapiens_assembly38.fasta"
NORMAL_BAM=$( ls /cluster/jxfu/proj/CIDC/Sentieon/data/mutation/${cohort}/*/normal_recalibrated.bam | while read s;do echo "-i ${s}";done)
PADDING=0
OUT_PON="/cluster/jxfu/proj/CIDC/Sentieon/data/PoN/${cohort}"
release_dir="/cluster/jxfu/proj/CIDC/Sentieon/release/sentieon-genomics-201808.01"
logfile="pipe_log/pon_${cohort}.log"
#0.Setup
export SENTIEON_LICENSE=172.24.216.24:8990
exec >$logfile 2>&1

#1. Create a panel of normal
echo "----------Create a panel of norma----------"
$release_dir/bin/sentieon driver -t $NUMBER_THREADS -r $REFERENCE \
$NORMAL_BAM --algo CNV \
--target $TARGET_BED \
--target_padding $PADDING --create_pon $OUT_PON
