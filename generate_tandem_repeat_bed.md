# Generate Tandem Repeat Annotation .bed File

## Overview
* pbsv Github page recommends including a tandem repeat annotation .bed \
during the discovery phase to "increase sensitivity and recall"
* I'm starting with trying Tandem Repeat Finder (trf)
  * I'm tried using a wrapper called trfbig that runs trf on large sequences, \
but it did not work properly
  * I installed both trf and trfbig (usc-trfbig) into by `NGS_analysis` \
environment via the bioconda channel
* Final File:
  * `/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes/poplar_14_5_v1_tandemrepeat.bed`
## Test run trying trfBig
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/tmp
trfBig /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta trf_test_1.fa
```
* The trfBig wrapper is not working properly

## Test out the regular trf program
* try with different size files

```
trf /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta 2 7 7 80 10 50 2000 -d

head -100000 /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta > poplar_14.5_top100klines.fasta

trf poplar_14.5_top100klines.fasta 2 7 7 80 10 50 2000 -h -d

sed -n '1,10000p; 10001q' /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta > poplar_14.5_Chr01_top10k_lines.fasta

sed -n '634531,644531p; 644532q' /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta > poplar_14.5_Chr02_top10k_lines.fasta

cat poplar_14.5_Chr01_top10k_lines.fasta \
poplar_14.5_Chr02_top10k_lines.fasta > poplar_14.5_Chrs_1_2_toplines.fasta

trf poplar_14.5_Chrs_1_2_toplines.fasta 2 7 7 80 10 50 2000 -h -d

grep -n '>' /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta

```

## Divide reference into chromosomes
* I think it will be faster and easier to process trf output from \
individual chromosomes
* Therefore, I need to divide the reference into chromosomes
* Will use `sed` to create the chromsome files and generate the `sed` \
commands in R
### Generate file with line number for start of each chromosome
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
grep -n '>' /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta > pop_14.5_chrom_start_lines.txt
wc -l /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/Populus_trichocarpa_var_14.5.mainGenome.fasta
```
  * total number of lines is 4903236
### Generate `sed` commands using R
```
chr_line_file <- paste('/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes', 'pop_14.5_chrom_start_lines.txt', sep = '/')
chr_lines <- read.table(chr_line_file, sep = '\t', header = F, 
  stringsAsFactors = F)
chr_lines <- rbind(chr_lines, '4903236:>last_line')
chr_split <- data.frame(matrix(
  data = unlist(strsplit(chr_lines[,1], split = ':')), ncol = 2, byrow = T),
  stringsAsFactors = F)
chr_split[,2] <- gsub('>', '', chr_split[,2])
chr_split[,3] <- as.numeric(chr_split[,1]) - 1
chr_split[nrow(chr_split),3] <- as.numeric(chr_split[nrow(chr_split), 1])
chr_split[,4] <- NA
chr_split[c(1:(nrow(chr_split)-1)),4] <- chr_split[c(2:nrow(chr_split)),3]
chr_split[,5] <- chr_split[,4] + 1
#
chr_split_sed_coms <- paste('sed -n \'', chr_split[,1], ',', chr_split[,4], 
  'p; ', chr_split[,5], 'q\' /home/f1p1/tmp/PBSV/Poplar14.5/REFERENCE/', 
  'Populus_trichocarpa_var_14.5.mainGenome.fasta > poplar_14.5_v1_', 
  chr_split[,2], '.fasta', sep = '')
split_coms_2 <- chr_split_sed_coms[-length(chr_split_sed_coms)]
split_coms_3 <- gsub('; 4903237q', '', split_coms_2)

chr_split_sed_com_file <- paste('/home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes', 'split_ref_into_chroms.sh', sep = '/')
write.table(split_coms_3, file = chr_split_sed_com_file, quote = F, sep = '\t',
  row.names = F, col.names = F)
```
### Make chromsome files for the reference
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
bash split_ref_into_chroms.sh
```

## run Tandem Repeat Finder (trf) on each chromosome
### Generate qsub files
* copied old qsub files to generate `trf_14.5_Chr01.sh` and \
`trf_14.5_scaffold_120.sh` files
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
bash
for i in {02..19};
do sed 's/01/'"$i"'/g' trf_14.5_Chr01.sh > trf_14.5_Chr$i.sh ; done

for i in 2190 2269 3504 3526 758;
do sed 's/120/'"$i"'/g' trf_14.5_scaffold_120.sh > trf_14.5_scaffold_$i.sh;
done

```
### submit jobs
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
bash
for i in {01..19};
do qsub trf_14.5_Chr$i.sh; done

for i in trf_14.5_scaffold_*.sh;
do qsub $i; done
```
### Backup output
```
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
mkdir trf_chrom_dat_files
cp *.dat ./trf_chrom_dat_files
tar -cvzf trf_chrom_dat.tar.gz ./trf_chrom_dat_files
```

## Process trf output to generate .bed
### Remove excess lines
```
bash
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
for i in {01..19};
do sed -e '1,15d' 'poplar_14.5_v1_Chr'$i'.fasta.2.7.7.80.10.50.2000.dat' > \
'Chr'$i'no_header.dat';
done

for i in 120 2190 2269 3504 3526 758;
do sed -e '1,15d' 'poplar_14.5_v1_scaffold_'$i'.fasta.2.7.7.80.10.50.2000.dat' \
> 'scaffold_'$i'no_header.dat';
done
```
### Extract columns 1,2,3, and 8
* 1 = start point
* 2 = end point
* 3 = size of repeat - will use to generate name
* 8 = aliignment score
```
bash
cd /home/f1p1/tmp/poplar_branches/ref_stuff/poplar_var_14.5_V1_chromosomes
for i in {01..19};
do cut -d " " -f1-3,8 'Chr'$i'no_header.dat' > 'Chr'$i'cols_1_2_3_8.dat';
done

for i in 120 2190 2269 3504 3526 758;
do cut -d " " -f1-3,8 'scaffold_'$i'no_header.dat' > \
'scaffold_'$i'cols_1_2_3_8.dat';
done
```
### Generate .bed file in R
* Requirements
  * Column 01 = chromosome name
  * Column 02 = repeat start position
  * Column 03 = repeat end position
  * Column 04 = name of repeat
    * My plan is to add 'mer' to the end of column 3 from `..1_2_3_8.dat` to \
to give a sense of the length of the repeat
  * Column 05 = score
    * My plan is to use Column 8 which is the allignment score
  * Column 06 = strand
    * My plan is to use '+' for all of them
#### Commands
```
data_dir <- paste('/home/f1p1/tmp/poplar_branches/ref_stuff/', 
  'poplar_var_14.5_V1_chromosomes/', sep = '')
col_files <- system(paste('ls ', data_dir, '*_1_2_3_8.dat', sep = ''), 
  intern = T)

for(cf in col_files){
  tmp_df <- read.table(cf, header = F, sep = ' ', 
    stringsAsFactors = F)
  tmp_df$rpt_name <- paste(tmp_df[,3], 'mer', sep = '')
  tmp_chr_name <- gsub('cols_1_2_3_8.dat', '', 
    gsub(paste('/home/f1p1/tmp/poplar_branches/ref_stuff/', 
    'poplar_var_14.5_V1_chromosomes/', sep = ''), '', cf))
  tmp_df$chr_name <- tmp_chr_name
  tmp_df$strand <- '+'
  tmp_bed <- data.frame(chrom = tmp_df$chr_name, start = tmp_df[, 1],
    end = tmp_df[, 2], rpt_name = tmp_df$rpt_name, score = tmp_df[, 4],
    stand = tmp_df$strand, stringsAsFactors = F)
  tmp_bed_file <- paste(data_dir, tmp_chr_name, '_tandemrepeat.bed.tmp', 
    sep = '')
  write.table(tmp_bed, file = tmp_bed_file, quote = F, sep = '\t', 
    row.names = F, col.names = F)
}
```
### Combine files into single .bed file
```
cat *.bed.tmp > poplar_14_5_v1_tandemrepeat.bed
```


* next: generate chr name files, repeat name files, and + files for columns\
to generate the .bed

* need to figure out how to process the output to generate a .bed
  * I think I want to remove the header and then extract the first few/
columns using a linux command, then generate some sort of name and score for
each element, and then convert it into a .bed file





  * Figure out how to generate a .bed file from the output of trf

* Chr02 starts on 634531
