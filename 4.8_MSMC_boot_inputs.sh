#!/bin/bash -l
#SBATCH --job-name=bootdata
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 01:00:00
#SBATCH --mem=8GB
#SBATCH -p bml
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_boot_inputs.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_boot_inputs.err

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset

INDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/inputs/phased

conda activate msmc2

# get list of chr files
INPUTS=$(find ${INDIR}/all10_*_msmc.txt | paste -sd " ")

mkdir -p $INDIR/boot
cd $INDIR/boot

# generates 200GB of sequence
~/bin/msmc-tools/multihetsep_bootstrap.py -n 20 --nr_chromosomes 24 bootstrap_number ${INPUTS}
