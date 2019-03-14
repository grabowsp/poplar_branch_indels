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
* I tried to incorporate separate `...subreads.bam` files into a single \
`pbmm2 align` command but couldn't get it to work
* It's not clear how to properly map the `...subreads.bam` files individually \
and then merge the outputs
* So will merge the `...subreads.bam` files into a single `.bam` file for \
mapping
* will use bamtools `merge` function for that
  * requires making a file of file names of the `...subreads.bam` files
### PAXL
* Just testing that `bamtools merge` will work
  * Take-home: it workds, but it takes a long time to run, so should submit \
the job on the cluster - see PAXN
```
bash
source activate NGS_analysis
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
ls m*subreads.bam > PAXL_subreads_files.fofn
bamtools merge -list PAXL_subreads_files.fofn -out PAXL.all.subreads.bam

```
### PAXN
* testing that running it via submision works properly
  * Take-home: it works
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXN
ls m*subreads.bam > PAXN_subreads_files.fofn
qsub PAXN.mergesubreads.sh
```
### Generate file of subread file names for each library
* already manually generated them for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK;
  do cd ./$LIB
  ls m*subreads.bam > $LIB'_subreads_files.fofn'
  cd ..; done
```
### Generate submit files
* already merged the files for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK;
do sed 's/PAXN/'"$LIB"'/g' ./PAXN/PAXN.mergesubreads.sh \
> ./$LIB/$LIB.mergesubreads.sh; done
```
### Submit merge jobs
* already merged the files for PAXL and PAXN
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS
for LIB in PAYK;
do cd ./$LIB
qsub $LIB.mergesubreads.sh
cd ..; done
```

## Align PacBio reads to refernce genome
### Overview
* use `pbmm2` to map `..subreads.bam` files to reference
  * had to merge the `subreads.bam` files to get this to work properly
* NOTE: in future, need to adjust script so output is written to tmp \
directory and transferred to target directory when done
### Example code
```
source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/NGS_analysis

pbmm2 align \
/home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta \
PAXL.all.subreads.bam ref.PAXL_v2_r1.bam --sort --sample 'PAXL' --median-filter
```
### PAXL
* Check that command words
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.align_v2_r1.sh
```
* Make submit scripts for reps 2 and 3
```
for i in _r2 _r3;
do 
sed 's/_r1/'"$i"'/g' PAXL.align_v2_r1.sh > PAXL.align_v2$i.sh;
done
```
* submit jobs
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.align_v2_r2.sh
qsub PAXL.align_v2_r3.sh
```
### PAXN and PAYK
* Generate submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.align_v2_r1.sh > \
  ./$LIB/$LIB.align_v2_r1.sh
  for i in _r2 _r3;
    do sed 's/_r1/'"$i"'/g' ./$LIB/$LIB.align_v2_r1.sh > \
    ./$LIB/$LIB.align_v2$i.sh;
  done;
done
```
* Submit jobs
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do cd ./$LIB
  for i in _r1 _r2 _r3;
#    do echo $LIB.align_v2$i.sh;
    do qsub $LIB.align_v2$i.sh;
  done;
  cd ..;
done
```

## Discover signatures of structural variation
### Overview
* use `pbsv discover` to generate `...svsig.gz` file
* GitHub page highly recommends using a tandem repeat annotation .bed file
* NOTE: in future, need to adjust code so writes output to tmp directory \
and then transfer output to target directory on completion
### Example code
```
pbsv discover --tandem-repeats \
/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed \
ref.PAXL_v2_r1.bam Ptr145v1.PAXL_v2_r1.svsig.gz
```
### PAXL
* make sure command works
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.disc_v2_r1.sh

```
* Make submit scripts for reps 2 and 3
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
for i in _r2 _r3;
do 
sed 's/_r1/'"$i"'/g' PAXL.disc_v2_r1.sh > PAXL.disc_v2$i.sh;
done
```
* submit jobs
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.disc_v2_r2.sh
qsub PAXL.disc_v2_r3.sh
```
### PAXN and PAYK
* Generate submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.disc_v2_r1.sh > \
  ./$LIB/$LIB.disc_v2_r1.sh
  for i in _r2 _r3;
    do sed 's/_r1/'"$i"'/g' ./$LIB/$LIB.disc_v2_r1.sh > \
    ./$LIB/$LIB.disc_v2$i.sh;
  done;
done
```
* Submit jobs
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do cd ./$LIB
  for i in _r1 _r2 _r3;
#    do echo $LIB.disc_v2$i.sh;
    do qsub $LIB.disc_v2$i.sh;
  done;
  cd ..;
done
```

## Call structural vairants and assign genotypes
### Overview
* use `pbsv call` to generate `...vcf` file
### Example Code
```
pbsv call /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta \
Ptr145v1.PAXL_v2_r1.svsig.gz Ptr145v1.PAXL_v2_r1.vcf
```
### PAXL
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.call_v2_r1.sh
```
* adjusted commands to deal with I/O issues from before
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
qsub PAXL.call_v2_r2.sh
```
* generate `..._r3.sh` file
```
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/PAXL
sed 's/_r2/_r3/g' ./PAXL.call_v2_r2.sh > ./PAXL.call_v2_r3.sh
qsub PAXL.call_v2_r2.sh
```
### PAXN and PAYK
* Generate submit scripts
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.call_v2_r2.sh > \
  ./$LIB/$LIB.call_v2_r2.sh
  for i in _r1 _r3;
    do sed 's/_r2/'"$i"'/g' ./$LIB/$LIB.call_v2_r2.sh > \
    ./$LIB/$LIB.call_v2$i.sh;
  done;
done
```
* Submit jobs
```
bash
cd /home/f1p1/tmp/PBSV/Poplar14.5/LIBS/
for LIB in PAXN PAYK;
  do cd ./$LIB
  for i in _r1 _r2 _r3;
#    do echo $LIB.disc_v2$i.sh;
    do qsub $LIB.call_v2$i.sh;
  done;
  cd ..;
done
```
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

