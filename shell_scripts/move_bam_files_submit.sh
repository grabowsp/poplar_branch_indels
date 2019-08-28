#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -V
#$ -cwd
#$ -N transfer_bams
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae

cd /home/f1p1/tmp/poplar_branches/lib_mapping/

/usr/bin/rsync -avuP ./*.bam /home/grabowsky_scratch/poplar_branch_files/

