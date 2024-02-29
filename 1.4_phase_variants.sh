#!/bin/bash -l
#SBATCH --job-name=phase
#SBATCH --array=1-24%15
#SBATCH --time 06:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/shapeit_phase_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/shapeit_phase_%a.err

# statistical phasing for use with MSMC2

# this runs per chromosome--very fast (~10 min/chr)
# which I then concatenated with bcftools concat

conda activate shapeit

BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/autosomes.bed
CHR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $BED | cut -f1)

cd /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34

shapeit4 --input bcf_variant_filterA_ac1.dp3.mis20_n34_masked.vcf.gz --region $CHR --output phased/filterA_${CHR}_phased.vcf.gz --sequencing
