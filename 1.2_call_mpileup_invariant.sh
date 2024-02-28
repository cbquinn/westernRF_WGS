#!/bin/bash -l
#SBATCH --job-name=mpileup
#SBATCH --array=1-27%15
#SBATCH --time 4-00:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/mpileup_invariant_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/mpileup_invariant_%a.err

# Call with only the 27 samples >10x. Instead of calling by chromosome, call by individual...
# awk '$2 >= 10 {print $0}' /group/ctbrowngrp2/cbquinn/fox4/0_bams/samplelist_34_wdepths.txt > /group/ctbrowngrp2/cbquinn/fox4/0_bams/samplelist_27_wdepths10xplus.txt

# Use these vcfs to calculate genome-wide heterozygosity in 1-mb windows

start=`date +%s`
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset

WORKDIR=/group/ctbrowngrp2/cbquinn/fox4

LIST=/group/ctbrowngrp2/cbquinn/fox4/0_bams/samplelist_27_wdepths10xplus.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f1)
MEAN_DP=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $LIST | cut -f2)
INDIR=/group/ctbrowngrp2/cbquinn/fox4/0_bams/indelrealigned
OUTDIR=/group/ctbrowngrp3/scratch/cbquinn/fox/vcfs/invariant
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed
REF=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/GCF_018345385.1_ASM1834538v1_genomic.fna

MAX_DP=$(($MEAN_DP * 2))


echo "Calling variant and invariant sites for $SAMPLE"

# only call for individuals in ms
# and only call for autosomes (with mappability mask applied)
# remove the -v flag to include invariant sites also

module load bcftools/1.14

bcftools mpileup -f $REF \
    --max-depth $MAX_DP \
    -q 30 --min-BQ 30 \
    -a DP,AD,ADF,ADR,INFO/AD \
    $INDIR/${SAMPLE}.bam | \
bcftools call -G - -m -Ob -o $OUTDIR/bcf_invariant_unfiltered_${SAMPLE}.bcf
bcftools index $OUTDIR/bcf_invariant_unfiltered_${SAMPLE}.bcf

# filter for minDP
bcftools filter \
	-i 'INFO/DP < 5' \
	-R $BED \
	-O b -o $OUTDIR/bcf_invariant_filtered_${SAMPLE}.bcf \ 
	$OUTDIR/bcf_invariant_unfiltered_${SAMPLE}.bcf
bcftools index $OUTDIR/bcf_invariant_filtered_${SAMPLE}.bcf

# remove unfiltered file (I don't think this worked...)
if [[ -s $OUTDIR/bcf_invariant_filtered_${SAMPLE}.bcf ]]
then
	echo "something went wrong..."
else
    rm $OUTDIR/bcf_invariant_unfiltered_${SAMPLE}.bcf*
fi

end=`date +%s`
runtime=$((end-start))
echo "finished!"
echo $runtime
