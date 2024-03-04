#!/bin/bash -l
#SBATCH --job-name=msmc
#SBATCH --array=1-10
#SBATCH --time 4-00:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/MSMC_vcf_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/MSMC_vcf_%a.err

STARTTIME=$(date +%s)

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset

# uses tab delimited list with 3 columns: sampleID, popID, meandepth:
# LAS_F01 LAS     28

LIST=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/msmc.bamlist
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f1)
POP=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f2)
DEPTH_MEAN=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f3)
DEPTH_MIN=$(($DEPTH_MEAN / 3))
DEPTH_MAX=$(($DEPTH_MEAN * 2))
REF=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/GCF_018345385.1_ASM1834538v1_genomic.fna
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/inputs
INDIR=/group/ctbrowngrp2/cbquinn/fox4/0_bams/indelrealigned

conda activate msmc2
module load htslib
module load samtools/1.14
module load bcftools/1.14

cd $OUTDIR

# outputs 2 files: an all-site vcf and a mask

for CHR in NC_054824.1 NC_054825.1 NC_054826.1 NC_054827.1 NC_054828.1 NC_054829.1 NC_054830.1 NC_054831.1 NC_054832.1 NC_054833.1 NC_054834.1 NC_054835.1 NC_054836.1 NC_054837.1 NC_054838.1 NC_054839.1 NC_054840.1 NC_054841.1 NC_054842.1 NC_054843.1 NC_054844.1 NC_054845.1 NC_054846.1 NC_054847.1

do

echo "creating vcf and mask for $SAMPLE at $CHR"

samtools mpileup -q 30 -Q 30 -C 50 -u -r $CHR -f $REF $INDIR/${SAMPLE}.bam | \
bcftools call -c -V indels | \
~/bin/msmc-tools/bamCaller.py $DEPTH_MEAN ${SAMPLE}_${CHR}_DEPTHmask.bed.gz | \
gzip -c > ${SAMPLE}_${CHR}_msmc.vcf.gz

done
