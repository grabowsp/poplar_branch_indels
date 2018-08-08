Steps used in analysis of indel variation in poplar branches


# Make VCFs that only contain chromosomes
When I first tried to compare the vcfs from the different samples, different
files would contain/miss indels on different non-anchored contigs, so needed
to make vcf files that only contained positions anchored on chromosomes

This is a shell script to use vcftools to generate the chromosome-only
.vcf files

`cd /home/grabowsky/tools/workflows/poplar_branch_indels/shell_scripts`
`bash make_chrom_only_vcfs.sh`


 
