#!/bin/bash -l
#SBATCH --job-name=boot
#SBATCH --array=1-20%3
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --time 6-00:00:00
#SBATCH --mem=80GB
#SBATCH -p bml
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_boot_run_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_boot_run_%a.err

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset


INDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/inputs/phased
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/outputs/phased

conda activate msmc2

cd $OUTDIR

p="10*1+22*2+1*4+1*6"
dir="t20_10x1_22x2_1x4_1x6"

# estimate Ne
# 64 time segments with 34 free parameters
mkdir -p $OUTDIR/$dir/boot

# n = 2
msmc2 -t 16 -i 20 -p $p -I 0,1,2,3 -o ${OUTDIR}/$dir/boot/LAS_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 4,5,6,7 -o ${OUTDIR}/$dir/boot/ORC_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 8,9,10,11 -o ${OUTDIR}/$dir/boot/RM_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 12,13,14,15 -o ${OUTDIR}/$dir/boot/WAC_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 16,17,18,19 -o ${OUTDIR}/$dir/boot/EAST_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
# n = 1
msmc2 -t 16 -i 20 -p $p -I 20,21 -o ${OUTDIR}/$dir/boot/RUS_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 22,23 -o ${OUTDIR}/$dir/boot/AK_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 24,25 -o ${OUTDIR}/$dir/boot/SN_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt
msmc2 -t 16 -i 20 -p $p -I 16,17 -o ${OUTDIR}/$dir/boot/ECAN_bootstrap_number_${SLURM_ARRAY_TASK_ID}.out ${INDIR}/boot/bootstrap_number_${SLURM_ARRAY_TASK_ID}/*txt


