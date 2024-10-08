#!/bin/bash -l
#SBATCH --job-name=filter
#SBATCH --time 24:00:00
#SBATCH --mem=1GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/filter_variant.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/filter_variant.err

# tested 4 different sets of filters and examined PCAs and bcftools/ROH
# opted to not use a maf because only one individual from AK and RUS, so distorts their heterozygosity and PCA/admixture relationships 
# (as a consequence, many of the 16 mil SNPs are Eurasia vs. North America--an only N. Amer set would be more like 12 mil

WORKDIR=/group/ctbrowngrp2/cbquinn/fox4
cd $WORKDIR

FILE=$WORKDIR/1_vcfs/variant/n41/unfiltered/bcf_variant_unfiltered_42_autosome.bcf
OUTDIR=$WORKDIR/1_vcfs/variant/n34
BEDmask=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed

module load bcftools/1.14
module load bedtools2/2.30.0

# FILTERS APPLIED
# remove known admixed individuals that are not part of this study
# max pooled depth filter of 1150 [sum(2*meanIndDp)]
# min pooled depth filter of 191 [sum(meanIndDp/3)]
# remove snps within 5 bp of indels
# site QUAL score >100 (this is very conservative, to make up for no maf)
# at least 6 reads to support the minor allele (across indivs)
# remove sites with excess het (>80% indiv) to get rid of paralogs 
# heterozygous genotypes to have at least 20% of the total reads supporting the minor allele and homozygous genotypes to have <10% of reads containing a different allele
# biallelic, remove indels
# <20% missigness

# Filters that vary (A, B, C, D)
# min depth to call genotype = 3 or 5
# allele count = 1 or 3 


#################
### FILTER A ####
#################

# AC=1; MINDP=3; MISS=0.2  

# change sample names
bcftools reheader -s /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n41/newnames_41.txt $FILE | \
# remove known admixed individuals not in study
bcftools view -S /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/samplenames_34.txt -Ou | \
# set genotypes with too few reads to missing
bcftools filter -S . -e 'FMT/DP< 3' -Ou | \
# set genotypes with allelic balance discrepancies to missing
bcftools filter -S . -e ' GT="HOM"  & (sMin(FORMAT/AD) / FORMAT/DP) > 0.2 ' | \
bcftools filter -S . -e ' GT="HET"  & (sMin(FORMAT/AD) / FORMAT/DP) < 0.2 ' | \
# remove sites with pooled depths too high or too low
bcftools view  -e 'INFO/DP > 1150 | INFO/DP <191' | \
# remove sites within 5 bp of indels
bcftools filter -g 5 -Ou | \
# remove sites that don't have at least 5 reads (across) all samples supporting the minor allele
bcftools view -i 'INFO/AD[1]>5 && INFO/AD[0]>5' | \
# remove sites that are not biallelic (more or less than 2 alleles) and have a site quality score < 30
bcftools view -m2 -M2 -v snps -c 1:minor -i '%QUAL>=100' -Ou | \
# remove sites where >80% of individuals are heterozygous (or >70% of individuals when >10% are missing)
bcftools filter -e ' (F_PASS(GT="het") >= 0.80) || ( F_PASS(GT="het") >= 0.70 && F_MISSING>0.1) ' -Ou | \
#bcftools filter -e 'F_PASS(GT="het") >= 0.80 || MQBZ < -3 ' -Ou | \
# remove sites with too many missing individuals
bcftools filter -e 'F_MISSING>0.2' -Oz \
-o $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_autosome.vcf.gz

# index
bcftools index -f $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_autosome.vcf.gz

# mask
bedtools intersect \
	-a $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_autosome.vcf.gz \
	-b $BEDmask \
	-wa -header | bcftools view -Oz -o $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked.vcf.gz 
bcftools index -f $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked.vcf.gz 

# get table of missigness per individual
bcftools stats -s - $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked.vcf.gz | grep -E ^PSC | cut -f3,14 > $OUTDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked.imiss

#################
### FILTER B ####
#################

# AC=1; MINDP=5; MISS=0.2

# change sample names
bcftools reheader -s /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n41/newnames_41.txt $FILE | \
# remove known admixed individuals not in study
bcftools view -S /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/samplenames_34.txt -Ou | \
# set genotypes with too few reads to missing
bcftools filter -S . -e 'FMT/DP< 5' -Ou | \
# set genotypes with allelic balance discrepancies to missing
bcftools filter -S . -e ' GT="HOM"  & (sMin(FORMAT/AD) / FORMAT/DP) > 0.2 ' | \
bcftools filter -S . -e ' GT="HET"  & (sMin(FORMAT/AD) / FORMAT/DP) < 0.2 ' | \
# remove sites with pooled depths too high or too low
bcftools view  -e 'INFO/DP > 1150 | INFO/DP <191' | \
# remove sites within 5 bp of indels
bcftools filter -g 5 -Ou | \
# remove sites that don't have at least 5 reads (across) all samples supporting the minor allele
bcftools view -i 'INFO/AD[1]>5 && INFO/AD[0]>5' | \
# remove sites that are not biallelic (more or less than 2 alleles) and have a site quality score < 30
bcftools view -m2 -M2 -v snps -c 1:minor -i '%QUAL>=100' -Ou | \
# remove sites where >80% of individuals are heterozygous (or >70% of individuals when >10% are missing)
bcftools filter -e ' (F_PASS(GT="het") >= 0.80) || ( F_PASS(GT="het") >= 0.70 && F_MISSING>0.1) ' -Ou | \
#bcftools filter -e 'F_PASS(GT="het") >= 0.80 || MQBZ < -3 ' -Ou | \
# remove sites with too many missing individuals
bcftools filter -e 'F_MISSING>0.2' -Oz \
-o $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_autosome.vcf.gz

# index
bcftools index -f $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_autosome.vcf.gz

# mask
bedtools intersect \
	-a $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_autosome.vcf.gz \
	-b $BEDmask \
	-wa -header | bcftools view -Oz -o $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_masked.vcf.gz 
bcftools index -f $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_masked.vcf.gz 

# get table of missigness per individual
bcftools stats -s - $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_masked.vcf.gz | grep -E ^PSC | cut -f3,14 > $OUTDIR/bcf_variant_filterB_ac1.dp5.mis20_n34_masked.imiss

#################
### FILTER C ####
#################

# AC=3; MINDP=3; MISS=0.2

# change sample names
bcftools reheader -s /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n41/newnames_41.txt $FILE | \
# remove known admixed individuals not in study
bcftools view -S /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/samplenames_34.txt -Ou | \
# set genotypes with too few reads to missing
bcftools filter -S . -e 'FMT/DP< 3' -Ou | \
# set genotypes with allelic balance discrepancies to missing
bcftools filter -S . -e ' GT="HOM"  & (sMin(FORMAT/AD) / FORMAT/DP) > 0.2 ' | \
bcftools filter -S . -e ' GT="HET"  & (sMin(FORMAT/AD) / FORMAT/DP) < 0.2 ' | \
# remove sites with pooled depths too high or too low
bcftools view  -e 'INFO/DP > 1150 | INFO/DP <191' | \
# remove sites within 5 bp of indels
bcftools filter -g 5 -Ou | \
# remove sites that don't have at least 5 reads (across) all samples supporting the minor allele
bcftools view -i 'INFO/AD[1]>3 && INFO/AD[0]>3' | \
# remove sites that are not biallelic (more or less than 2 alleles) and have a site quality score < 30
bcftools view -m2 -M2 -v snps -c 3:minor -i '%QUAL>=100' -Ou | \
# remove sites where >80% of individuals are heterozygous (or >70% of individuals when >10% are missing)
bcftools filter -e ' (F_PASS(GT="het") >= 0.80) || ( F_PASS(GT="het") >= 0.70 && F_MISSING>0.1) ' -Ou | \
#bcftools filter -e 'F_PASS(GT="het") >= 0.80 || MQBZ < -3 ' -Ou | \
# remove sites with too many missing individuals
bcftools filter -e 'F_MISSING>0.2' -Oz \
-o $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_autosome.vcf.gz

# index
bcftools index -f $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_autosome.vcf.gz

# mask
bedtools intersect \
	-a $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_autosome.vcf.gz \
	-b $BEDmask \
	-wa -header | bcftools view -Oz -o $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_masked.vcf.gz 
bcftools index -f $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_masked.vcf.gz 

# get table of missigness per individual
bcftools stats -s - $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_masked.vcf.gz | grep -E ^PSC | cut -f3,14 > $OUTDIR/bcf_variant_filterC_ac3.dp3.mis20_n34_masked.imiss

#################
### FILTER D ####
#################

# AC=3; MINDP=5; MISS=0.2

# change sample names
bcftools reheader -s /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n41/newnames_41.txt $FILE | \
# remove known admixed individuals not in study
bcftools view -S /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/samplenames_34.txt -Ou | \
# set genotypes with too few reads to missing
bcftools filter -S . -e 'FMT/DP< 5' -Ou | \
# set genotypes with allelic balance discrepancies to missing
bcftools filter -S . -e ' GT="HOM"  & (sMin(FORMAT/AD) / FORMAT/DP) > 0.2 ' | \
bcftools filter -S . -e ' GT="HET"  & (sMin(FORMAT/AD) / FORMAT/DP) < 0.2 ' | \
# remove sites with pooled depths too high or too low
bcftools view  -e 'INFO/DP > 1150 | INFO/DP <191' | \
# remove sites within 5 bp of indels
bcftools filter -g 5 -Ou | \
# remove sites that don't have at least 5 reads (across) all samples supporting the minor allele
bcftools view -i 'INFO/AD[1]>3 && INFO/AD[0]>3' | \
# remove sites that are not biallelic (more or less than 2 alleles) and have a site quality score < 30
bcftools view -m2 -M2 -v snps -c 3:minor -i '%QUAL>=100' -Ou | \
# remove sites where >80% of individuals are heterozygous (or >70% of individuals when >10% are missing)
bcftools filter -e ' (F_PASS(GT="het") >= 0.80) || ( F_PASS(GT="het") >= 0.70 && F_MISSING>0.1) ' -Ou | \
#bcftools filter -e 'F_PASS(GT="het") >= 0.80 || MQBZ < -3 ' -Ou | \
# remove sites with too many missing individuals
bcftools filter -e 'F_MISSING>0.2' -Oz \
-o $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_autosome.vcf.gz

# index
bcftools index -f $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_autosome.vcf.gz

# mask
bedtools intersect \
	-a $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_autosome.vcf.gz \
	-b $BEDmask \
	-wa -header | bcftools view -Oz -o $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_masked.vcf.gz 
bcftools index -f $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_masked.vcf.gz

# get table of missigness per individual
bcftools stats -s - $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_masked.vcf.gz | grep -E ^PSC | cut -f3,14 > $OUTDIR/bcf_variant_filterD_ac3.dp5.mis20_n34_masked.imiss
```


### tabulate number of variants
and also look over missingness per individual 

```
paste *imiss | cut -f1,2,4,6,8 | sort -k 2 > imiss_filters_masked.txt

rm nsites_filters_autosome.txt
for i in $(ls *autosome.vcf.gz)
do
SAMPLE=$(echo $i)
VARIANTS=$(bcftools index -n $i)
echo -e "$SAMPLE\t$VARIANTS" >> nsites_filters_autosome.txt
done

rm nsites_filters_masked_genmap.txt
for i in $(ls *masked.vcf.gz)
do
SAMPLE=$(echo $i)
VARIANTS=$(bcftools index -n $i)
echo -e "$SAMPLE\t$VARIANTS" >> nsites_filters_masked_genmap.txt
done

# note: these are the final counts for the different filters 
# filterID | filter | n_autosome | n_mask_map | n_mask_mapRP
# --- | --- | --- | --- | ---
# A | ac1.dp3.mis20 | 16,922,841 | 16,654,053
# B | ac1.dp5.mis20 | 16,163,355 | 15,977,905
# C | ac3.dp3.mis20 | 12,620,284 | 12,415,459
# D | ac3.dp5.mis20 | 11,826,897 | 11,692,506
# the mappability filter on its own removes about 200K snps

