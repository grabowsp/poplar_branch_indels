#!/bin/bash

# Script for calling the different R scripts for generating the stats for
#   pbsv outputs

source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_Indel_analysis.r $1
Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_BND_analysis.r $1
Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_DUP_analysis.r $1
Rscript /home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/pbsv_INV_analysis.r $1


