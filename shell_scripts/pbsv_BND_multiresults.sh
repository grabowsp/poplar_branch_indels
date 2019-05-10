#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -V
#$ -cwd
#$ -N pbsv_BND_multisamp

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_ngmlr_r03.r

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_ngmlr_r04.r

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_pbmm2_r01.r

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_pbmm2_r02.r

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_pbmm2_r03.r

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_v2.2_combo_pbmm2_r04.r

