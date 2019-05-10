#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -V
#$ -cwd
#$ -N pbsv_BND_ngmlr_r2

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_ngmlr_r02.r


