#!/bin/bash -l
#SBATCH --job-name=trim
#SBATCH --array=1-73%25
#SBATCH -p med
#SBATCH -A ctbrowngrp
#SBATCH --mem 1GB 
#SBATCH -t 2-00:00:00
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/split_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/split_%a.err

### script uses a tab-delimited LIST file that  formatted like this:
#JOB	RUN	RFID	SID	BASENAME	LIB	PLATFORM
#1	run1	RF1	S11-0008	RF1_S11-0008_run1	lib1_SI	hiseq
#2	run1	RF2	S15-0575	RF2_S15-0575_run1	lib1_SI	hiseq

### trim reads

# Some tracking info
start=`date +%s`
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Define variables
LIST=/group/ctbrowngrp3/scratch/cbquinn/fox/fastqlist.txt
INDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/fastqs/untrimmed/
OUTDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/fastqs/trimmed/

# check col numbers
BASENAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f5)

echo "Running trimgalore on ${INDIR}${BASENAME}"

module load fastqc
module load cutadapt
module load trim_galore/1

trim_galore --gzip -o $OUTDIR --paired ${INDIR}${BASENAME}_R1.fastq.gz ${INDIR}${BASENAME}_R2.fastq.gz --length 30

end=`date +%s`
runtime=$((end-start))
echo "finished trimming:"
echo $runtime

### split reads into smaller chunks for parallel alignment

echo "Splitting trimmed fastq for read $BASENAME into chunks of 100-million read"
mkdir -p ${OUTDIR}split
zcat ${OUTDIR}${BASENAME}*_1.fq.gz | split -a 4 -l 400000000 -d --additional-suffix=".fastq" - ${OUTDIR}split/${BASENAME}.R1
zcat ${OUTDIR}${BASENAME}*_2.fq.gz  | split -a 4 -l 400000000 -d --additional-suffix=".fastq" - ${OUTDIR}split/${BASENAME}.R2

gzip ${OUTDIR}split/${BASENAME}*fastq

end=`date +%s`
runtime=$((end-start))
echo "finished splitting:"
echo $runtime
