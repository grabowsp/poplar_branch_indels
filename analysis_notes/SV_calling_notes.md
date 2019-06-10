# Detecting and calling SV in poplar data using poplar ref v4.0
## Overview
### Goals
* Call SVs in branches via multiple methods
* Identify method(s) that call SVs with enough confidence to map accumulation
of somatic SV mutations going up the tree
* Validate some SVs using Illumina data
* Describe "what are the SV mutations?" (general goal of project)
* Characterize where are the mutations found and do they affect genes? (general
goal of project)
### Steps
1) Call SVs using HAGSC "Old" PacBio pipeline
  * each sample separate (default) and all samples together
2) Call SVs using HAGSC "New" PacBio pipeline
  * each sample separate (default) and all samples together
3) Call SVs using ngmlr/SNIFFLES/SURVIVOR
  * NOTE: as of 10/17/2018, there seems to be a bug in SURVIVOR with collecting read-support \
data which is causing issues with aggregating individual VCFs, calling genotypes, etc.
  * need to wait for update/correction before re-run pipeline and look at genotypes in combined \
files
4) Call SVs using VarScan
5) General/summary statistics
5) Generate statistics that indicate the quality of SV calling
6) Compare the results from each method
7) Generate pipeline for checking Illumina data
8) Validate 10 SVs in Illunina data (as test case)

## PacBio "Old" pipeline
* Mike automatically ran this pipeline
* Notes about transfering files to $WORK directory: `./pacbio_SVcaller_notes.md`
### VCF files
* Found in `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/old_PB_SVcaller/`
* Individual Library files: `pop_branches_LIBNAME__old_PB_SVcaller_Ref14.5_v1.0.vcf`
  * ex: `pop_branches_PAXL_old_PB_SVcaller_Ref14.5_v1.0.vcf`
* Combined .vcf: `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/old_PB_SVcaller/pop_branches_combo_old_PB_SVcaller_Ref14.5_v1.0.vcf`

## PacBio "New" pipeline
* Chris ran this pipeline for me
* Notes about transfering files to $WORK directory: `./pacbio_SVcaller_notes.md`
### VCF files
* Found in `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/`
* Individual Library files: `ref.LIBNAME.vcf`
  * ex: `ref.PAXL.vcf`
* Combined .vcf: `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/ref.ALLData.vcf`

## ngmlr/SNIFFLES/SURVIVOR
### Steps
1) ngmlr on each library separately using "SAM attributes enables"
2) sort .bams
3) SNIFFLES on each .bam separately
4) sort vcfs
5) SURVIVOR on all vcfs
6) rerun SNIFFLES to force-call all SVs based on SURVIVOR output
7) sort vcfs
8) Re-run SURVIVOR to get merged vcf
### Analysis Notes
`./ngmlr_sniffles_notes.md`
### VCF files
* Found in `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/sniffles`
* Individual "raw" genotypes calls: `LIBNAME.14.5v1.0Ref.ngmlr.sorted.sniffles.sorted.vcf`
  * ex: `PAXL.14.5v1.0Ref.ngmlr.sorted.sniffles.sorted.vcf`
* Combined .vcf: `/home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/sniffles/popbranch.ngmlr.sniffles.survivor.supervisedmerged1kbdist.vcf`


