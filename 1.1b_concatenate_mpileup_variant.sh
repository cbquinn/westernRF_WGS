#!/bin/bash -l
#SBATCH --job-name=concat
#SBATCH --time 02:00:00
#SBATCH --mem=500M
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/concat_variant.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/concat_variant.err

# Concatenate bychr files into one vcf

start=`date +%s`
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset

WORKDIR=/group/ctbrowngrp2/cbquinn/fox4

cd $WORKDIR/1_vcfs/variant/n41/unfiltered/

module load bcftools/1.14

# indexes vcf files
FILENAMES=$(ls bcf_variant_unfiltered_42_chr{1..24}.bcf)
for i in $FILENAMES
do
    bcftools index $i
 done

bcftools concat -o bcf_variant_unfiltered_42_autosome.bcf -O b ${FILENAMES}
bcftools index bcf_variant_unfiltered_42_autosome.bcf

end=`date +%s`
runtime=$((end-start))
echo "finished concatenating:"
echo $runtime
