cd /home/t4c1/WORK/grabowsk/data/poplar_branches/SV_calling_analysis/new_PB_SVcaller/chrom_only

for i in `ls *.vcf`; do vcf-sort $i | bgzip > $i'.gz'; done

for i in `ls *.gz`; do tabix -p vcf $i; done
