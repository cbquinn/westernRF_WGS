#!/bin/bash -l
#SBATCH --job-name=het
#SBATCH --array=1-27%15
#SBATCH --time 2-00:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/window_het_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/window_het_%a.err

LIST=/group/ctbrowngrp2/cbquinn/fox4/0_bams/samplelist_27_wdepths10xplus.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f1)
INDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/vcfs/invariant
BCF=$INDIR/bcf_invariant_filtered_${SAMPLE}.bcf
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/autosomes_1MBwindows.bed
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/2_heterozygosity/windows/1MB/

module load bcftools/1.14
module load bedtools2/2.30.0

cd $OUTDIR

for CHR in NC_054824.1 NC_054825.1 NC_054826.1 NC_054827.1 NC_054828.1 NC_054829.1 NC_054830.1 NC_054831.1 NC_054832.1 NC_054833.1 NC_054834.1 NC_054835.1 NC_054836.1 NC_054837.1 NC_054838.1 NC_054839.1 NC_054840.1 NC_054841.1 NC_054842.1 NC_054843.1 NC_054844.1 NC_054845.1 NC_054846.1 NC_054847.1

do

sed -n "/^$CHR/p" $BED > ${SAMPLE}.${CHR}.bed

# create bed file of hets and then intersect with windows bedfile
bcftools view -r $CHR $BCF | \
bcftools query -f '%CHROM\t%POS\t[%GT]\n' | \
awk '$3 == "0/1" {print}' | \
awk -v OFS="\t" '{print $1, $2-1, $2, $3}' | \
bedtools intersect -a ${SAMPLE}.${CHR}.bed -b - -c -sorted > ${SAMPLE}.${CHR}.hets.txt

# create bed file of non-missing sites and then intersect with windows bedfile
# very large "-b" files cause bedtools to soak up huge amounts of memory. Try running this by itself with increased memory and with the "sorted" algorithm
bcftools view -r $CHR $BCF | \
bcftools query -f '%CHROM\t%POS\t[%GT]\n' | \
awk '$3 != "./." {print}' | \
awk -v OFS="\t" '{print $1, $2-1, $2, $3}' | \
bedtools intersect -a ${SAMPLE}.${CHR}.bed -b - -c -sorted > ${SAMPLE}.${CHR}.sites.txt 

paste ${SAMPLE}.${CHR}.hets.txt ${SAMPLE}.${CHR}.sites.txt | awk -v OFS="\t" '{print $1, $2, $3, $4, $8}' > ${SAMPLE}.${CHR}.combined.tsv

done

cat ${SAMPLE}.*.combined.tsv > ${SAMPLE}.genomewidehet_windows.tsv
