#!/bin/bash -l
#SBATCH --job-name=bwa_mod
#SBATCH --array=1-89%25
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 3-00:00:00
#SBATCH --mem=8GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/bwamem_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/bwamem_%a.err


# make things fail on errors
set -o nounset
set -o errexit

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# LIST is a text file with as many rows as there are split fastq chunkks and the following columns: JOB READ1 READ2 CHUNK BASENAME RUN RFID SID LIB PLATFORM POP
# use these to add detailed readgroup info

LIST=/group/ctbrowngrp3/scratch/cbquinn/fox/fastqs/trimmed/split/splitfastqlist.txt
READ1=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f2)
READ2=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f3)
CHUNK=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f4)
BASENAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f5)
RUN=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f6)
RFID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f7)
SID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f8)
LIB=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f9)
PLATFORM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f10)
POP=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f11)

REF=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/GCF_018345385.1_ASM1834538v1_genomic.fna
INDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/fastqs/trimmed/split/
OUTDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/bams/

echo "Aligning paired chunk ${CHUNK} to arctic fox reference using bwa-mem"

module load bwa
module load samtools

# Align R1 and R2 
# Pipe to a bam file that excludes bad mapping scores
# Pipe to a sorted bam for merging

cd $OUTDIR
bwa mem ${REF} ${INDIR}${READ1} ${INDIR}${READ2} | \
samtools view -h -b - | \
samtools sort -o ${OUTDIR}chunk/${CHUNK}.bam -

echo "adding readgroups..."
module load picardtools/2.7.1

java -jar ${PICARD}/picard.jar AddOrReplaceReadGroups \
    I=${OUTDIR}chunk/${CHUNK}.bam \
    O=${OUTDIR}chunk/readgroups/${CHUNK}.bam \
    RGID=${RUN} RGLB=${LIB} RGPL=${PLATFORM} RGPU=NULL RGSM=${RFID}_${SID} \
    VALIDATION_STRINGENCY=LENIENT

# remove raw output bam
rm ${OUTDIR}chunk/${CHUNK}.bam
### NOTE! 
fastqs for ORC_S18-2071_run4 & WAC_S10-0583_run4 got swapped. I corrected at the bam stage by renaming. The fastqs uploaded to NCBI are correct.

