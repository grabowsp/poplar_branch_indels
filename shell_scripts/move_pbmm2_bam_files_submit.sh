#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -V
#$ -cwd
#$ -N transfer_pbmm2_bams
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae

cd /home/f1p2/tmp/Poplar14.5_pbsv/pbsv_v2.1.1_run

/usr/bin/rsync -avuP ./*.bam /home/grabowsky_scratch/poplar_branch_files/pbmm2_bams

