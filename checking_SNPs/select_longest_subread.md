# Notes about selecting the longest subread from the alignments

## Overview
* Want to select the longest subread (pass from the same sequencing molecule) \
so that same read doesn't over-contribute to genotype designations

## Strategy
* sort bam by readname
  * `samtools sort -n`
* Extract readnames using 'view'
  * might need to break up into chromosomes
* Extract read length using 'length' (?)
* Choose reads to keep
  * `view -R FILE` where FILE has RG (read-group) tags to keep - not sure \
this will work

## Test
* Need to find the sort script I used to get this to work
```
cd /home/grabowsky_scratch/poplar_branch_files

```
* adjust this script to sort a bam by readname
```
#!/bin/bash

#$ -q "all.q"
#$ -V
#$ -pe smp 8
#$ -cwd
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -N PBAU_sortbam

source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

samtools sort -m 2G -o PBAU.14.5v1.0Ref.ngmlr.sorted.bam -@ 8 PBAU.14.5v1.0Ref.ngmlr.bam


```




