#!/bin/bash -l
#SBATCH --job-name=fsc_ORCLAS
#SBATCH --time 3-00:00:00
#SBATCH --mem=2GB
#SBATCH --array=1-10
#SBATCH --cpus-per-task=2
#SBATCH -A ctbrowngrp
#SBATCH -p bmm
#SBATCH -o /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240625/slurmlogs/fastsimcoal_sim_revision_RMORC_%a.out
#SBATCH -e /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240625/slurmlogs/fastsimcoal_sim_revision_RMORC_%a.err

fsc=~/bin/fsc27_linux64/fsc27093
cd /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240625/out

echo "run ${SLURM_ARRAY_TASK_ID}: for ${1}"

cp ../in/${1}.tpl ../in/${1}_${SLURM_ARRAY_TASK_ID}.tpl
cp ../in/${1}.est ../in/${1}_${SLURM_ARRAY_TASK_ID}.est

$fsc -t ../in/${1}_${SLURM_ARRAY_TASK_ID}.tpl -n 1 -e ../in/${1}_${SLURM_ARRAY_TASK_ID}.est -E10000 -s 10000 -d -q -I -B2 -c2

