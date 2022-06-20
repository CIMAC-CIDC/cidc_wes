#!/bin/bash
# Len Taing 2022 (TGBTG)
# script to merge consensus CNV regions of a particular type (either GAIN or LOSS)
# takes two args- arg $1 is the file, arg $2 is the type {GAIN, LOSS}

#First check if there are regions of that type in the file
#For example, if we have all GAINS, then merging LOSS regions should just
#be an empty file
if grep -q $2 $1
then
    grep $2 $1 | mergeBed -c 4 -o distinct
fi

