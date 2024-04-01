#!/bin/bash -l
#SBATCH --job-name=fsc
#SBATCH --time 3-00:00:00
#SBATCH --mem=20GB
#SBATCH --array=1-30
#SBATCH --cpus-per-task=2
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/fastsimcoal_sim_%A_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/fastsimcoal_sim_%A_%a.err

# script takes argument $1 that is name of model (iso, iso.mig, iso.mig.bot)

fsc=~/bin/fsc27_linux64/fsc27093
cd /group/ctbrowngrp3/scratch/cbquinn/fox/abc/out_2

echo "run ${SLURM_ARRAY_TASK_ID}: for ${1}"

cp ../in_2/${1}.tpl ../in_2/${1}_${SLURM_ARRAY_TASK_ID}.tpl
cp ../in_2/${1}.est ../in_2/${1}_${SLURM_ARRAY_TASK_ID}.est

$fsc -t ../in_2/${1}_${SLURM_ARRAY_TASK_ID}.tpl -n 1 -e ../in_2/${1}_${SLURM_ARRAY_TASK_ID}.est -E10000 -s 10000 -d -q -I
