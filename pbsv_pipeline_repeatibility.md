# Test of repeatibility of the pbsv pipeline for detecting SVs in the poplar \
branches

## Overview
* Run recommended pbsv pipeline, as layed out on the pbsv Github page on \
Jan. 29 2019
* Run 3 times on 3 samples to look at repeatibility
  * PAXL
  * PAXN
  * PAYK
## Generate Tandem Repeat Annotation .bed
### Overview
* Github page recommends including a tandem repeat annotation .bed during \
the discovery phase to "increase sensitivity and recall"
### Info about generating .bed file
* `./generate_tandem_repeat_bed.md`
### .bed file
* `/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed`

## Discover signatures of structural variation
### Overview
* use `pbmm2` to map `..subreads.bam` files to reference
### Example code
```
source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

pbmm2 align /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta m*.subreads.bam m*.subreads.bam ref.PAXL_v2_r1.bam --sort --sample 'PAXL' --median-filter
```
### PAXL
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.align_v2_r1.sh
qsub PAXL.align_v2_r2.sh
qsub PAXL.align_v2_r3.sh
```

