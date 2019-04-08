# Rerun pbsv pipeline using v2.1.1 of pbsv (Feb 2019)

## Overview
* Run recommended pbsv v2.1.1 pipeline, as layed out on the pbsv Github page \
on Jan. 29 2019
* PAXL, PAXN, and PAYK had already been run as part of the repeatibility \
analysis

## Generate Tandem Repeat Annotation .bed
### Overview
* Github page recommends including a tandem repeat annotation .bed during \
the discovery phase to "increase sensitivity and recall"
* This info is also incuded in the .md file for repeatibility analysis
### Info about generating .bed file
* `./generate_tandem_repeat_bed.md`
### .bed file
* `/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed`

## Merge .subreads.bam files to make single file for mapping
### Overview
* I tried to incorporate separate `...subreads.bam` files into a single \
`pbmm2 align` command but couldn't get it to work
* It's not clear how to properly map the `...subreads.bam` files individually \
and then merge the outputs
* So will merge the `...subreads.bam` files into a single `.bam` file for \
mapping
* will use bamtools `merge` function for that
  * requires making a file of file names of the `...subreads.bam` files
### Generate file of subread file names for each library
* already manually generated files them for PAXL, PAXN and PAYK
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
ls m*subreads.bam > $LIB'_subreads_files.fofn'
cd ..; done
```
* In future, should use `ls $PWD/m*subreads.bam` so get full path of the files

### Generate submit files
* already merged the files for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXN/'"$LIB"'/g' ./PAXN/PAXN.mergesubreads.sh \
> ./$LIB/$LIB.mergesubreads.sh; done
```
### Submit merge jobs
* already merged the files for PAXL, PAXN and PAYK
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
qsub $LIB.mergesubreads.sh
cd ..; done
```

## Align PacBio reads to refernce genome
### Overview
* use `pbmm2` to map `..subreads.bam` files to reference
  * had to merge the `subreads.bam` files to get this to work properly
### Example code
```
source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

pbmm2 align \
/home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta \
PAXL.all.subreads.bam ref.PAXL_v2_r1.bam --sort --sample 'PAXL' --median-filter
```
### Rest of the Libraries
* Generate submit scripts
  * Already generated scripts for PAXL, PAXN, and PAYK
  * NOTE: for future runs, need to adjust script so that writes to tmp file \
and final files are transferred to the desired directory
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.align_v2_r1.sh > \
  ./$LIB/$LIB.align_v2_r1.sh;
done
```
* Submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
  do cd ./$LIB
  qsub $LIB.align_v2_r1.sh
  cd ..;
done
```
## Discover signatures of structural variation
### Overview
* use `pbsv discover` to generate `...svsig.gz` file
* GitHub page highly recommends using a tandem repeat annotation .bed file
### Example code
```
pbsv discover --tandem-repeats \
/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed \
ref.PAXL_v2_r1.bam Ptr145v1.PAXL_v2_r1.svsig.gz
```
### Generate submit scripts
* Already generated scripts for PAXL, PAXN, and PAYK
* Note: as doing this, adjusted script to deal with I/O overload
  * PAYZ and PAZF used old version of script used for PAXL, PAXN, and PAYK
  * rest of libraries used newer version of script
#### Scripts for PAYZ and PAZF
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.disc_v2_r1.sh > \
  ./$LIB/$LIB.disc_v2_r1.sh;
done
```
#### Scripts for remaining libraries
* adjust submit script to write to tmp directory then transfer to desired \
directory
* Testing for PAZG - works
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAZG
qsub PAZG.disc_v2_r1.sh
```
* change submit scripts for remaining libraries
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAZH PBAT PBAU PBAW;
  do sed 's/PAZG/'"$LIB"'/g' ./PAZG/PAZG.disc_v2_r1.sh > \
  ./$LIB/$LIB.disc_v2_r1.sh;
done
```
### Submit jobs
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
  do cd ./$LIB
  qsub $LIB.disc_v2_r1.sh
  cd ..;
done
```

## Call structural vairants and assign genotypes
### Overview
* use `pbsv call` to generate `...vcf` file
* NOTE: Next time, need to include the -j flag to set the number of threads
### Example Code
```
pbsv call /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta \
Ptr145v1.PAXL_v2_r1.svsig.gz Ptr145v1.PAXL_v2_r1.vcf
```
### Generate Submission Scripts
#### PAZY and PAZF
* these libraries were run with the previous approach so I want their \
submit files to accurately represent the previous commands
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF;
  do sed 's/PAYK/'"$LIB"'/g' ./PAYK/PAYK.call_v2_r1.sh > \
  ./$LIB/$LIB.call_v2_r1.sh;
done
```
#### Remaining libraries
* Manually generated file for PAZG
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAZH PBAT PBAU PBAW;
  do sed 's/PAZG/'"$LIB"'/g' ./PAZG/PAZG.call_v2_r1.sh > \
  ./$LIB/$LIB.call_v2_r1.sh;
done
```
### Submit jobs
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
  do cd ./$LIB
  qsub $LIB.call_v2_r1.sh
  cd ..;
done
```

## All Samples
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
qsub full.call_v2_try1.sh
```


* NEXT Steps:
  * DONE: check on all the mapping results
  * DONE: Check on PAXL SV calling/genotyping results
  * DONE: Generate commands for discovery for rest of the libraries
  * Submit discovery jobs
  * Check on discovery results
  * Generate commands for sv calling for individuals
  * Submit sv calling for individuals jobs
  * Generate commands fof sv calling using all samples
    * Generate 4 with different sample orders
  * Submit sv calling for all samples jobs

## Compare new individual results to previous pipeline
### Overview
* First will `compare ref.PAXL.vcf` (old pipeline) to \
`Ptr145v1.PAXL_v2_r1.vcf` (new pipeline)
* Things I'm interested:
  * Overlap in SVs
    * `vcftools --diff-site`
  * Differences in coverage/depth
    * I don't get why the aligned .bam files were so much bigger for the \
second round - perhaps Chris can help explain that
  * Relative numbers of each SV class
  * Size distribution in each SV class


## Move big .bam files to new directory but set up symlinks
### try something like this
* `PAXL_big_bam_names.txt`
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
NEW_DIR=/home/f1p2/tmp/Poplar14.5_pbsv/pbsv_v2.1.1_run
for BAMFILE in `ls ref.PAXL_v2*bam*`;
  do echo ln -s $NEW_DIR/$BAMFILE $BAMFILE >> PAXL_symlink_coms.txt; done
/usr/bin/rsync -avuP ref.PAXL_v2*bam* $NEW_DIR
rm ref.PAXL_v2*bam*
bash PAXL_symlink_coms.txt
```
### Submit job test for PAXN
* run after the 'discovery' script has run
* `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXN/transfer_run2_bams_PAXN.sh`
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXN
qsub transfer_run2_bams_PAXN.sh
```
* worked
### Make scripts for rest of libraries
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXN/'"$LIB"'/g' ./PAXN/transfer_run2_bams_PAXN.sh > \
./$LIB/transfer_run2_bams_$LIB.sh;
done
```
### Submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
qsub transfer_run2_bams_$LIB.sh
cd ..;
done
```
