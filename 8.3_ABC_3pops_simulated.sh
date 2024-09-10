#!/bin/bash -l
#SBATCH --job-name=fsc_meta2
#SBATCH --time 3-00:00:00
#SBATCH --mem=5GB
#SBATCH --array 1-50%30
#SBATCH --cpus-per-task=2
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240615meta/slurmlogs/fastsimcoal_sim_revision_meta2_%a.out
#SBATCH -e /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240615meta/slurmlogs/fastsimcoal_sim_revision_meta2_%a.err

fsc=~/bin/fsc27_linux64/fsc27093
cd /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240719meta/out

echo "run ${SLURM_ARRAY_TASK_ID}: for ${1}"

cp ../in/${1}.tpl ../in/${1}_${SLURM_ARRAY_TASK_ID}.tpl
cp ../in/${1}.est ../in/${1}_${SLURM_ARRAY_TASK_ID}.est

$fsc -t ../in/${1}_${SLURM_ARRAY_TASK_ID}.tpl -n 1 -e ../in/${1}_${SLURM_ARRAY_TASK_ID}.est -E 1000 -s 10000 -d -q -I -c2 -B2

