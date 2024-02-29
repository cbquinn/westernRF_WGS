#!/bin/bash -l
#SBATCH --job-name=split_psmcfa
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 00:10:00
#SBATCH --mem=1GB
#SBATCH -p high
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/split_psmcfa.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/split_psmcfa.err

start=`date +%s`

for i in 1 2 5 12
do

LIST=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/PSMC/pscm.bamlist
SAMPLE=$(sed "${i}q;d" $LIST | cut -f2)
INDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/PSMC

cd $INDIR

echo "splitting psmc.fa files for ${SAMPLE} for bootstrapping"

/group/ctbrowngrp2/hennelly/hennelly/bin/psmc/utils/splitfa fas/${SAMPLE}.psmc.fa > fas/split/${SAMPLE}_split.psmc.fa


done

end=`date +%s`
runtime=$((end-start))
echo "finished splitting psmc.fa for ${SAMPLE} bootstrapping after ${runtime} seconds"
