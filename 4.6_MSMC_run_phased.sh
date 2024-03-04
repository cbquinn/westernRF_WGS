#!/bin/bash -l
#SBATCH --job-name=MSMC
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --time 3-00:00:00
#SBATCH --mem=60GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_run.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_run.err

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset
set -e

INDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/inputs/phased
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/outputs/phased

conda activate msmc2

# get list of chr files
INPUTS=$(find $INDIR/all10_*_msmc.txt | paste -sd " ")

# 64 time segments with 34 free parameters
p="10*1+22*2+1*4+1*6"
dir="t20_10x1_22x2_1x4_1x6"

# "The `-s` flag tells MSMC to skip sites with ambiguous phasing. As a rule of thumb: For population size estimates, we have found that unphased sites are not so much of a problem, but for cross-population analysis we typically remove those." -- https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md

# estimate Ne
mkdir -p $OUTDIR/$dir
msmc2 -t 16 -i 20 -p $p -I 0,1,2,3 -o ${OUTDIR}/$dir/LAS_msmc.out $INPUTS
msmc2 -t 16 -i 20 -p $p -I 4,5,6,7 -o ${OUTDIR}/$dir/ORC_msmc.out $INPUTS
msmc2 -t 16 -i 20 -p $p -I 8,9,10,11 -o ${OUTDIR}/$dir/RM_msmc.out $INPUTS
msmc2 -t 16 -i 20 -p $p -I 12,13,14,15 -o ${OUTDIR}/$dir/WAC_msmc.out $INPUTS
msmc2 -t 16 -i 20 -p $p -I 16,17,18,19 -o ${OUTDIR}/$dir/EAST_msmc.out $INPUTS


# Print out final statistics about resource use before job exits
scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
