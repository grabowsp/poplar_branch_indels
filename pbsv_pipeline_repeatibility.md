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

## Merge .subreads.bam files to make single file for mapping
### Overview
* I'm not sure how to map the ...subreads.bam files individually, so will /
merge the ...subreads.bam files into a single .bam file for mapping
* will use bamtools `merge` function for that
### PAXL
```
bash
source activate NGS_analysis
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
ls m*subreads.bam > PAXL_subreads_files.fofn
bamtools merge -list PAXL_subreads_files.fofn -out PAXL.all.subreads.bam

```
### PAXN
```
bash
ls m*subreads.bam > PAXN_subreads_files.fofn
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXN
qsub PAXN.mergesubreads.sh
```

### Generate file of subread file names for each library
* already manually generated them for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do ls ./$LIB/m*subreads.bam > $LIB_subreads_files.fofn; done  
```
### Generate submit files
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXN/'"$LIB"'/g' ./PAXN/PAXN.mergesubreads.sh \
> ./$LIB/$LIB.mergesubreads.sh; done
```
* NEXT: wait to make sure A) PAXN merge worked, and B) PAXL alignment words/\
makes a difference, then submit these jobs to merge bams

## Discover signatures of structural variation
### Overview
* use `pbmm2` to map `..subreads.bam` files to reference
* NOTE: I need to re-do this - I need to combine the subread files and then/
map them
### Example code
```
source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

pbmm2 align /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta m*.subreads.bam ref.PAXL_v2_r1.bam --sort --sample 'PAXL' --median-filter
```
### PAXL
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.align_v2_r1.sh
qsub PAXL.align_v2_r2.sh
qsub PAXL.align_v2_r3.sh
```

## Discover signatures of structural variation
### Overview
* use `pbsv discover` to generate `...svsig.gz` file
### Example code

### PAXL

