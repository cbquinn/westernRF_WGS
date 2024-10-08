
#!/bin/bash -l
#SBATCH --job-name=mlrho
#SBATCH --array=1-34
#SBATCH --time 24:00:00
#SBATCH --mem=500M
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/mlrho_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/mlrho_%a.err

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounsetq

module load gsl

WORKDIR=/group/ctbrowngrp2/cbquinn/fox4

# I run this using an array. With a tab-delimited file ("$LIST") with <sampleID> <mean depth> for each individual. Depth is an integer.
# Also uses a bed file ("$BED") of regions in autosomes (excludes any unplaced scaffolds, X chr, and mitochondrion) and "good" mappability (e.g., excludes repetitive regions or regions with low mapping score)
# Should take ~12 hours for a 30x genome, uses very little memory

LIST=/group/ctbrowngrp2/cbquinn/fox4/0_bams/samplelist_34_wdepths.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f1)
MEAN_DP=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f2)
INDIR="$WORKDIR/0_bams/indelrealigned"
OUTDIR="$WORKDIR/2_heterozygosity/mlrho"
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed

MAX_DP=$(($MEAN_DP * 2))

# Options for mlrho:
# -c <NUM> minimum coverage; default: 4. Can change this depending on average coverage of samples
# -M maximum distance analyzed in rho computation; at 0 distance, all positions are looked at individually

module load samtools/1.14
cd $OUTDIR

# /home/hennelly/bin/MlRho_2.9

samtools view -L $BED -b $INDIR/$SAMPLE.bam | \
samtools mpileup -d $MAX_DP -q 30 -Q 30 - | \
cut -f 2,5 | \
awk -f ~/bin/MlRho_2.9/Scripts/bam2pro.awk | \
~/bin/FormatPro_0.5/formatPro -c 5 -n ${SAMPLE}.c5
# compute the genome-wide mutation rate
/home/hennelly/bin/MlRho_2.9/mlRho -M 0 -n ${SAMPLE}.c5 > theta.${SAMPLE}.c5.txt
