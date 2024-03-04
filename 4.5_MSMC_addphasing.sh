#!/bin/bash -l
#SBATCH --job-name=MSMCinputPH
#SBATCH --array=1-10
#SBATCH --nodes 1
#SBATCH --time 01:00:00
#SBATCH --mem=4GB
#SBATCH -p bml
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_addphasing_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_addphasing_%a.err

LIST=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/msmc.bamlist
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f1)

INDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/inputs
phasedIN=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/phased/bfilterA_phased.vcf.gz

cd $INDIR

# Inputs: 
# (1) phased and filtered vcf for all individuals
# (2) bam2vcf file that has already been made for this chr

module load bcftools/1.14

# make temp file
echo $SAMPLE > ${SAMPLE}.name

for CHR in NC_054824.1 NC_054825.1 NC_054826.1 NC_054827.1 NC_054828.1 NC_054829.1 NC_054830.1 NC_054831.1 NC_054832.1 NC_054833.1 NC_054834.1 NC_054835.1 NC_054836.1 NC_054837.1 NC_054838.1 NC_054839.1 NC_054840.1 NC_054841.1 NC_054842.1 NC_054843.1 NC_054844.1 NC_054845.1 NC_054846.1 NC_054847.1
do


# split the phased vcf by individual and chr
# Pipe to: Extract sites in the unphased vcf, not present in the phased sites (i.e. those sites which could not be phased)

bcftools view -s $SAMPLE -r $CHR -O z $phasedIN > temp/${SAMPLE}_${CHR}_phased.vcf.gz
tabix temp/${SAMPLE}_${CHR}_phased.vcf.gz


# bgzip and reheader and index msmcIN
gunzip -c ${SAMPLE}_${CHR}_msmc.vcf.gz | bgzip | \
bcftools reheader -s ${SAMPLE}.name | bcftools view -o temp/${SAMPLE}_${CHR}_msmc.vcf.gz -Oz
tabix temp/${SAMPLE}_${CHR}_msmc.vcf.gz

# get sites that were not able to be phased so that they can be added
bcftools isec -n -1 -p temp/${SAMPLE}_${CHR} -c all temp/${SAMPLE}_${CHR}_msmc.vcf.gz temp/${SAMPLE}_${CHR}_phased.vcf.gz

# Concatenate unphased sites with phased sites and sort the vcf.
# 0 is sites private to vcf A (unphased sites)
bcftools concat temp/${SAMPLE}_${CHR}_phased.vcf.gz temp/${SAMPLE}_${CHR}/0000.vcf  | \
bcftools sort -O z -o phased/${SAMPLE}_${CHR}_phased_final.vcf.gz

# get sites that were not able to be phased so that they can be added
bcftools isec -n -1 -p temp/${SAMPLE}_${CHR} -c all temp/${SAMPLE}_${CHR}_msmc.vcf.gz temp/${SAMPLE}_${CHR}_phased.vcf.gz

# Concatenate unphased sites with phased sites and sort the vcf.
# 0 is sites private to vcf A (unphased sites)
bcftools concat temp/${SAMPLE}_${CHR}_phased.vcf.gz temp/${SAMPLE}_${CHR}/0000.vcf  | \
bcftools sort -O z -o phased/${SAMPLE}_${CHR}_phased_final.vcf.gz

# ^^^ phased and unphased sites for an indiv/chr

# remove temp files
rm -r temp/${SAMPLE}_${CHR}
rm temp/${SAMPLE}_${CHR}*

done
