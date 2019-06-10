# Run Poplar Branch PacBio data through pbsv v2.2

## Overview
* Updated pbsv from v2.1.1 to v2.2 on March 7, 2019
* Will run pbsv using:
  * alignment using pbmm2 - the PacBio aligner
    * ran previously for the v2.1.1 pipeline
    * For this alignment, I used the subread .bam files rather than the \
subread .fasta.gz files that I used for ngmlr. So, it won't exactly be a \
apple-to-apple comparison. I suppose I could re-run pbmm2 using the fasta.gz \
files, but it seems like pbmm2 is designed for the subread.bam files, so \
using those files should produce the best quality results...in theory
    * Lori had issues with alignments using this alogrithm for human data
  * alignment using ngmlr
    * Ran previously for the sniffles/SURVIVOR attempts
### Strategy
* Generate new parent directory
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs`
* Generate new subdirectories for each branch
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for LIB in PAXL PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW; do
  mkdir $LIB; done
```
* Save output using pbmm2 and ngmlr outputs to same directory
* use following for naming
  * include 'pbmm2' or 'ngmlr' to signify which alignement was used
  * use 'PtStettler14' to signify that branch 14.5 reference was used

## Add Read Group info to ngmlr alignment .bam files
* tried running `pbsv discover` on the output from ngmlr and got an error \
saying: `ERROR: BamHeader: read group ID not found`
* In future runs of ngmlr, I think the `--rg-id` flag will help fix that
* For now, I will add the the read group info afterwards
* Using picard to add read group info
  * using the `AddOrReplaceReadGroups` program in picard to add the read group \
info
### Test interactive session with PAXL
* takes about 1hr to run - should run remotely for remaining libraries
```
picard AddOrReplaceReadGroups \
I=/home/f1p1/tmp/poplar_branches/lib_mapping/PAXL.14.5v1.0Ref.ngmlr.sorted.bam \
O=/home/f1p1/tmp/poplar_branches/lib_mapping/PAXL.14.5v1.0Ref.ngmlr.sorted.withRG.bam \
RGID=PAXL_combined \
RGLB=PAXL \
RGPL=PACBIO \
RGPU=default_unit1 \
RGSM=PAXL
```
### Testing Submit scripts test with PAXN
* shell script
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXN/addngmlrReadGroups.PAXN.sh`
#### Submit job
```
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXN
qsub addngmlrReadGroups.PAXN.sh
```
* there will eventually be a new synatx, but not yet
### Rest of Libraries
#### Generate submit scripts
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXN/'"$LIB"'/g' ./PAXN/addngmlrReadGroups.PAXN.sh \
> ./$LIB/addngmlrReadGroups.$LIB.sh;
done
```
#### Submit jobs
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/
for LIB in PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd $LIB
qsub addngmlrReadGroups.$LIB.sh
cd ..;
done
```

## Discover Signatures of Structural Variation
### Test with PAXL:
* using pbmm2 alignement:
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PAXL.pbmm2.v2.2.disc.r01.sh`
* using ngmlr alignment:
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PAXL.ngmlr.v2.2.disc.r01.sh`
#### submit scripts
```
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL
qsub PAXL.pbmm2.v2.2.disc.r01.sh
qsub PAXL.ngmlr.v2.2.disc.r01.sh
```
### Rest of Libraries
#### Generate submit scripts
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for LIB in PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.pbmm2.v2.2.disc.r01.sh \
> ./$LIB/$LIB.pbmm2.v2.2.disc.r01.sh
sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.ngmlr.v2.2.disc.r01.sh \
> ./$LIB/$LIB.ngmlr.v2.2.disc.r01.sh;
done
```
#### Submit scripts
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for LIB in PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
qsub $LIB.pbmm2.v2.2.disc.r01.sh
qsub $LIB.ngmlr.v2.2.disc.r01.sh
cd ..;
done
```

## Call SVs
### Test with PAXL
* using pbmm2 alignment
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PAXL.pbmm2.v2.2.call.r01.sh`
* using ngmlr alignment
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL/PAXL.ngmlr.v2.2.call.r01.sh`
* submit scripts
```
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/PAXL
qsub PAXL.pbmm2.v2.2.call.r01.sh
qsub PAXL.ngmlr.v2.2.call.r01.sh
```
### Remaining libraries
#### Generate Scripts
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for LIB in PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.pbmm2.v2.2.call.r01.sh \
> ./$LIB/$LIB.pbmm2.v2.2.call.r01.sh
sed 's/PAXL/'"$LIB"'/g' ./PAXL/PAXL.ngmlr.v2.2.call.r01.sh \
> ./$LIB/$LIB.ngmlr.v2.2.call.r01.sh;
done
```
#### Submit jobs
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for LIB in PAXN PAYK PAYZ PAZF PAZG PAZH PBAT PBAU PBAW;
do cd ./$LIB
qsub $LIB.pbmm2.v2.2.call.r01.sh
qsub $LIB.ngmlr.v2.2.call.r01.sh
cd ..;
done
```

## Call genotypes in all libraries at once
### Replicate 1
* alphabetical order of libraries
* `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/full.call.pbmm2.v2.2.call.r01.sh`
* * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/full.call.ngmlr.v2.2.call.r01.sh`
```
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
qsub full.call.ngmlr.v2.2.call.r01.sh
qsub full.call.pbmm2.v2.2.call.r01.sh
```
### Rest of replicates
* rep 2, 3, and 4 ready to go
* using same sample orders as seen in `/home/f1p1/tmp/PBSV/Poplar14.5/LIBS`
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for SUBFILE in `ls *r02.sh`;
do qsub $SUBFILE; done
for SUBFILE in `ls *r03.sh`;
do qsub $SUBFILE; done
for SUBFILE in `ls *r04.sh`;
do qsub $SUBFILE; done
```

## Call genotypes in only 8 good branches
### Overview
* Will use same sample orders as before, but removing the two "bad" branches
* This will eventually give the most appropriate results
* Will omit PAZG (13.4) and PAYZ (14.1)
### Submision scripts
* example
  * `/home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs/full.call.ngmlr.v2.2.call.r01.8branch.sh`
### Submit
```
bash
cd /home/f1p1/tmp/poplar_branches/pbsv_v2.2_runs
for SUBFILE in `ls *8branch.sh`;
do qsub $SUBFILE; done
```

