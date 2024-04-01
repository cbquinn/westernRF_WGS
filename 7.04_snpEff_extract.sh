#!/bin/bash -l
#SBATCH --job-name=snpeff
#SBATCH --time 24:00:00
#SBATCH --mem=8GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/extract_snpeff.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/extract_snpeff.err

cd /group/ctbrowngrp2/cbquinn/fox4/5_load/snpeff

VCF=filterA_all_wGERP_AArotated_snpeff_all.vcf.gz

conda activate snpeff
module load bcftools
module load bedtools2

# get site list for each category
# remove sites with warnings

# create beds (with gerp scores) for each category

# intergenic - subset
SnpSift filter "ANN[0].EFFECT has 'intergenic_region' && !(exists ANN[0].ERRORS)" $VCF | \
bcftools query  -f '%CHROM\t%POS\n' | \
awk -v OFS="\t" '{print $1,$2-1, $2, "intergenic"}' |\
awk '{printf("%f\t%s\n",rand(),$0);}' | sort -t $'\t'  -T . -k1,1g | head -n 100000 | cut -f 2- | sort -k1,1 -k2,2 > intergenic.bed

SnpSift filter "ANN[0].EFFECT has 'intergenic_region' && !(exists ANN[0].ERRORS)" $VCF | \
bcftools view -Oz -o intergenic_full.vcf.gz

# synonymous
SnpSift filter "(ANN[0].EFFECT has 'synonymous_variant') && !(exists ANN[0].ERRORS)" $VCF | \
bcftools query  -f '%CHROM\t%POS\t%GERP\n' | \
awk -v OFS="\t" '{print $1,$2-1, $2, $3, "synonymous"}' > synonymous.bed

# missense
SnpSift filter "(ANN[0].EFFECT has 'missense_variant') && !(exists ANN[0].ERRORS)" $VCF | \
bcftools query  -f '%CHROM\t%POS\t%GERP\n' | \
awk -v OFS="\t" '{print $1,$2-1, $2, $3, "missense"}' > missense.bed

# will add in SIFT filtering here...

# lof
SnpSift filter "(exists LOF[0].PERC) && !(exists ANN[0].ERRORS)" $VCF | \
bcftools query  -f '%CHROM\t%POS\t%GERP\n' | \
awk -v OFS="\t" '{print $1,$2-1, $2, $3, "lof"}' > lof.bed

# missense
cat synonymous.bed missense.bed lof.bed > CDS.variants.bed

# subset to sites that have GERP scores
awk -F "\t" '{ if($4 != ".") { print $0} }' > CDS.variants.wGERP.bed

# create a vcf for each CDS category

for i in intergenic_full intergenic synonymous missense lof
do
echo "creating vcf for $i variants"
# create vcf
BED=${i}.bed
bedtools intersect -header \
    -a $VCF -b <(sortBed -i $BED) -wa | \
bcftools view -O z -o ${i}.vcf.gz
bcftools index ${i}.vcf.gz
done
