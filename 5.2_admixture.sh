#!/bin/bash -l
#SBATCH --job-name=admix
#SBATCH --array=1-15
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=4
#SBATCH --time 24:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/admixture_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/admixture_%a.err

# ALL SAMPLES
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/admixture/final/all
data="$INDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked_LDpruned.bed"

# REDUCED SAMPLES
#OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/admixture/final/reduced
#data="$INDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked_LDpruned_reduced.bed"

INDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/final
#program path
admix=/home/sophiepq/admixture_linux-1.3.0/admixture
K=$SLURM_ARRAY_TASK_ID
out=$(basename ${data%.bed})
module load deprecated/plink/1.90

cd $OUTDIR
mkdir -p K${K}
cd K${K}

echo "run   seed" > K${K}.random.txt

#runs 10 replicates
for run in {1..10};
do
    mkdir -p run${run}
    cd run${run}

    #set random seed
    seed=$RANDOM

    $admix --cv=10 $data $K -j4 -s $seed | tee run${run}.K${K}.out

    cd ..
    echo "run${run} $seed" >> K${K}.random.txt
done

scontrol show job ${SLURM_JOB_ID}
sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
