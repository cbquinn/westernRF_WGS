#!/bin/bash -l
#SBATCH --job-name=SIFT
#SBATCH --time 1-00:00:00
#SBATCH --mem=20GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/SIFT_annotate.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/SIFT_annotate.err

# Full instructions here: https://sift.bii.a-star.edu.sg/sift4g/AnnotateVariants.html

SIFT=/group/ctbrowngrp2/cbquinn/fox4/5_load/sift/SIFT4G_Annotator.jar
DATABASE=/group/ctbrowngrp3/scratch/cbquinn/fox/siftdb/Vulpes_lagopus/GCF_018345385
VCF=/group/ctbrowngrp2/cbquinn/fox4/5_load/snpeff/filterA_all_wGERP_AArotated_snpeff_all.vcf.gz
NAME=$(basename $VCF)
LOG=${NAME%.vcf.gz}_SIFT.log

conda activate sift
module load jdk/17.0.1

cd /group/ctbrowngrp2/cbquinn/fox4/5_load/sift/annotation_uniprot

gzip -cd ${VCF} > ${VCF%.gz}

java -Xmx16g -jar ${SIFT} -c -t -r ./ -d ${DATABASE} -i ${VCF%.gz} >> $LOG

#bgzip ${LOG%.log}$.vcf

# separate missense into tolerated and 
# prioritze most deleterious variant for multiple transcipts
# format as bedfile
module load bedtools2
module load bcftools 

grep "NONSYNONYMOUS" filterA_all_wGERP_AArotated_snpeff_all_SIFTannotations.xls | \
grep "DELETERIOUS" | \
sort -u -k1,1 -k2,2 | \
awk -v OFS="\t" '{print $1,$2-1, $2, $13, $14,  "missense-del"}' > SIFT_uniprot_missense_del.bed 
# 15770

grep "NONSYNONYMOUS" filterA_all_wGERP_AArotated_snpeff_all_SIFTannotations.xls | \
grep "TOLERATED" | \
sort -u -k1,1 -k2,2 | \
awk -v OFS="\t" '{print $1,$2-1, $2, $13, $14, "missense-tol"}' | \
# get records in a that are not in b
bedtools intersect -v -wa -a - -b SIFT_uniprot_missense_del.bed > SIFT_uniprot_missense_tol.bed 
# 26,123
# 41,893 total

# create a vcf for each CDS category

for i in missense_tol missense_del
do
echo "creating vcf for $i variants"
# create vcf
BED=SIFT_uniprot_${i}.bed
bedtools intersect -header \
    -a $VCF -b <(sortBed -i $BED) -wa | \
bcftools view -O z -o ${i}.vcf.gz
bcftools index ${i}.vcf.gz
done
