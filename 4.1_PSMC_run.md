#!/bin/bash -l
#SBATCH --job-name=psmc
#SBATCH --array=1-4
#SBATCH --time 3-00:00:00
#SBATCH --mem=3GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/PSMC_run_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/PSMC_run_%a.err

STARTTIME=$(date +%s)

# LIST is tab-delimited input file headed like this (only last 2 columns necessary):
# indelrealigned	RUS_S12-0237	13


LIST=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/PSMC/pscm.bamlist
BAMDIR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f1)
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f2)
DEPTH_MEAN=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f3)
DEPTH_MIN=$(($DEPTH_MEAN / 3))
DEPTH_MAX=$(($DEPTH_MEAN * 2))
REF=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/GCF_018345385.1_ASM1834538v1_genomic.fna
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/PSMC
INDIR=/group/ctbrowngrp2/cbquinn/fox4/0_bams

module load samtools/1.14
module load bcftools/1.14

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset


##### step 1: bam to fastq #####
# 2 days for bam2fq (<1 G), several hours for next part

# keep only autosomes with high mappability
# call all sites with mpileup
# convert from vcf to fasta

echo "creating fq for sample $SAMPLE and outputting file to $OUTDIR"

samtools view -L $BED -u $INDIR/$BAMDIR/${SAMPLE}.bam | \
samtools mpileup -Q 30 -q 30 -u -v -f $REF - | \
bcftools call -c | \
vcfutils.pl vcf2fq -d 5 -D $DEPTH_MAX -Q 30 | \
gzip > $OUTDIR/fqs/${SAMPLE}.fq.gz

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "bam to fastq done. Took $timestamp hours:minutes:seconds to complete..."

##### step 2: fastq to fasta #####

STARTTIME=$(date +%s)

/home/sophiepq/bin/psmc/utils/fq2psmcfa $OUTDIR/fqs/${SAMPLE}.fq.gz > $OUTDIR/fas/${SAMPLE}.psmc.fa

##### step 3: execute psmc #####

# test a few different parameter sets

t=15
p="4+25*2+4+6"
dir="t15_4_25x2_4_6"
mkdir -p $OUTDIR/outfiles/$dir
/home/sophiepq/bin/psmc/psmc -t $t -r5 -p $p -o  $OUTDIR/outfiles/$dir/${SAMPLE}.psmc $OUTDIR/fas/${SAMPLE}.psmc.fa

t=15
p="10*1+22*2+4+6"
dir="t15_10x1_25x2_4x2"
mkdir -p $OUTDIR/outfiles/$dir
/home/sophiepq/bin/psmc/psmc -t $t -r5 -p $p -o  $OUTDIR/outfiles/$dir/${SAMPLE}.psmc $OUTDIR/fas/${SAMPLE}.psmc.fa

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "psmc done. Took $timestamp hours:minutes:seconds to complete..."
