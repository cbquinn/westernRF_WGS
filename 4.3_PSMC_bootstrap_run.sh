#!/bin/bash -l
#SBATCH --job-name=split_psmcfa
#SBATCH --array=1-100%30
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 2-00:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/boot_psmc_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/boot_psmc_%a.err


ROUND=${SLURM_ARRAY_TASK_ID}
module load gnuplot

for SAMPLE in RUS_S12-0237 AK_S12-1159 ECAN_S14-0358 SN_S11-0008
do
STARTTIME=$(date +"%s")
echo "running run ${ROUND} of PSMC bootstrapping on ${SAMPLE}"

cd /group/ctbrowngrp2/cbquinn/fox4/4_SMC/PSMC

OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/PSMC
t=15
p="10*1+22*2+4+6"
dir="t15_10x1_25x2_4x2/boot/"
mkdir -p $OUTDIR/outfiles/$dir

seq 1 | xargs -i echo /home/sophiepq/bin/psmc/psmc -t $t -r5 -b -p $p -o $OUTDIR/outfiles/$dir/${SAMPLE}_boot${ROUND}.psmc $OUTDIR/fas/split/${SAMPLE}_split.psmc.fa | sh

cd $OUTDIR/outfiles/$dir/
/group/ctbrowngrp2/hennelly/hennelly/bin/psmc/utils/psmc_plot.pl -u 4.5e-09 -g 2 -R ${SAMPLE}_boot${ROUND} ${SAMPLE}_boot${ROUND}.psmc

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60,TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Took $timestamp hours:minutes:seconds to complete..."
done
