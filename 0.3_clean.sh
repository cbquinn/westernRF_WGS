#!/bin/bash -l
#SBATCH --job-name=clean
#SBATCH --array=1-34
#SBATCH --time 23:40:00
#SBATCH --mem=16GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/clean_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/clean_%a.err

STARTTIME=$(date +%s)

# make things fail on errors
set -o nounset
set -o errexit

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# LIST IS A TAB-DELIMITED LIST FORMATED WITH THE FOLLOWING HEADERS
# RFID    sampleID        SID     POP


LIST=/group/ctbrowngrp3/scratch/cbquinn/fox/bams/samplelist.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f2)
TEMPDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/bams/
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/0_bams/
REF=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/GCF_018345385.1_ASM1834538v1_genomic.fna
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/autosomes.bed

module load samtools/1.16.1
module load picardtools/2.7.1

cd $TEMPDIR

# merge split bams

INPUTS_LIST=$( cat <(ls ${TEMPDIR}chunk/readgroups/${SAMPLE}*.bam))
INPUTS=$(echo $INPUTS_LIST | tr "\n" " ")
INPUTS_NUMB=$(echo $INPUTS_LIST | wc -l)

echo "Samtools will now merge $INPUTS_NUMB files:"
echo "samtools merge -f ${TEMPDIR}merged/${SAMPLE}.bam $INPUTS"

samtools merge -f ${TEMPDIR}merged/${SAMPLE}.bam $INPUTS

# mark duplicates

echo "marking duplicates..."
mkdir -p deduped
java -jar ${PICARD}/picard.jar MarkDuplicates \
    VALIDATION_STRINGENCY=LENIENT \
      I=${TEMPDIR}merged/${SAMPLE}.bam \
      O=${TEMPDIR}deduped/${SAMPLE}.bam \
      M=${TEMPDIR}deduped/qcmetrics/${SAMPLE}_duplicationMETRICS.txt

# indel realignment

echo "indexing..."
samtools index ${TEMPDIR}deduped/${SAMPLE}.bam

module load GATK/3.6

java -Xmx4g -jar ${GATK_HOME}/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-nt 4 \
-R ${REF} \
-I ${TEMPDIR}deduped/${SAMPLE}.bam \
-o ${TEMPDIR}deduped/intervals/${SAMPLE}.intervals \
-log ${TEMPDIR}deduped/intervals/${SAMPLE}.intervals.log

java -Xmx4g -jar ${GATK_HOME}/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ${REF} \
-I ${TEMPDIR}deduped/${SAMPLE}.bam \
-targetIntervals ${TEMPDIR}deduped/intervals/${SAMPLE}.intervals \
-o ${OUTDIR}indelrealigned/${SAMPLE}.bam

echo "indexing..."
samtools index ${OUTDIR}indelrealigned/${SAMPLE}.bam

echo "getting stats from GATK..."
# collect allignment metrics
java -jar ${PICARD}/picard.jar CollectAlignmentSummaryMetrics \
          R=${REF} \
          I=${OUTDIR}indelrealigned/${SAMPLE}.bam \
          O=${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_summaryMETRICS.txt \
          VALIDATION_STRINGENCY=LENIENT

echo "getting stats from samtools..."

# count total reads
samtools view -L $BED -c ${OUTDIR}indelrealigned/${SAMPLE}.bam > ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_NoReads_tot.flagstat

# count number of good reads: exclude umapped reads, secondary and chimeric reads, singleton reads, and clones
samtools view -F 1796 -L $BED -c ${OUTDIR}indelrealigned/${SAMPLE}.bam > ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_NoReads_1796.flagstat

# count number of good reads with mapq > 30: exclude umapped reads, secondary and chimeric reads, singleton reads, and clones
samtools view -F 1796 -q 30 -L $BED -c ${OUTDIR}indelrealigned/${SAMPLE}.bam > ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_NoReads_1796q30.flagstat

# count number of reads that map to more than one location and/or have an umapped mate
samtools view -f 2312 -L $BED -c ${OUTDIR}indelrealigned/${SAMPLE}.bam > ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_NoReads_2312.flagstat

# count number of unmapped reads
samtools view -f 4 -L $BED -c ${OUTDIR}indelrealigned/${SAMPLE}.bam > ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_NoReads_4.flagstat

# echo "getting  histogram of mapping quality..."
# for "good" reads (co1 is mapQ, col2 is no. reads)
samtools view -F 1796 -L $BED ${OUTDIR}indelrealigned/${SAMPLE}.bam | cut -f 5 | sort | uniq  -c | sort -n | awk '{printf("%s\t%d\n",$2,$1);}' > ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_mapq.dist

echo "getting mean depth of coverage..."
# (default is to exclude reads with flag 1796)
conda activate mosdepth
mosdepth -n --fast-mode --by 50000 ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE} ${OUTDIR}indelrealigned/${SAMPLE}.bam
mosdepth -n --fast-mode --by 50000 --mapq 30 ${OUTDIR}indelrealigned/qcmetrics/${SAMPLE}_mapq30 ${OUTDIR}indelrealigned/${SAMPLE}.bam

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "All done. Took $timestamp hours:minutes:seconds to complete..."

