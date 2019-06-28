# Analysis of regions with large het SVs in Tree14

## Plan
1. Find SV's that don't have exact same name but are the same SV
1. Collect names of SV's greater than 50kb
  * Generate script for outputting filtered SV's above a certain size for \
each rep
  * Can repeat with smaller size classes later, if need be
2. Verify that SV's are real using IGV
  * 3 dels and 3 dups
2. Check individual VCF files to see if there is any indication about the \
quality of the SVs
3. Quantify the amount of genes in each SV
4. Quantify the amount of coding seq in each SV
5. Quantify the amount of intergenic seq in each SV
6. Calculate genome-wide patterns of gene content
  * Types:
    * gene density
    * coding sequence density
    * intergenic sequence density
  * Scales:
    * Genome-wide
    * Chromosome-wide
    * Within windows: 100kb, 50kb, 20kb
    * Non-centromeric and non-telomeric regions selected by window patterns

## Get names of Large SV's
### SV's > 50kb
* script used for getting SNP names
  * `/home/grabowsky/tools/workflows/poplar_branch_indels/r_scripts/large_SV_get_names.r`
      

