#!/bin/bash -l
#SBATCH --job-name=CCR
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --time 8-00:00:00
#SBATCH --mem=160GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_run.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_run.err

# fail on errors and don't allow unset variables
set –o pipefail
set –o errexit
set –o nounset
set -e

INDIR=/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/msmc2/inputs/phased
OUTDIR=/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/msmc2/outputs/phased

conda activate msmc2

# get list of chr files
INPUTS=$(find $INDIR/all13_*_msmc.txt | paste -sd " ")

p="10*1+22*2+1*4+1*6"

dir="t20_10x1_22x2_1x4_1x6"

# estimate cross-coalescence
# the -s flag tells MSMC to skip sites with ambiguous phasing. unphased sites are not so much of a problem for Ne but for cross-population analysis, it is recommended to remove them

#### n = 2


# LAS_EAST ccr
msmc2 -t 16 -i 20 -p $p -s -I 0-16,0-17,0-18,0-19,1-16,1-17,1-18,1-19,2-16,2-17,2-18,2-19,3-16,3-17,3-18,3-19 -o ${OUTDIR}/$dir/LAS_EAST_msmc.out $INPUTS

# ORC_EAST ccr
msmc2 -t 16 -i 20 -p $p -s -I 4-16,4-17,4-18,4-19,5-16,5-17,5-18,5-19,6-16,6-17,6-18,6-19,7-16,7-17,7-18,7-19 -o ${OUTDIR}/$dir/ORC_EAST_cross_msmc.out $INPUTS

# RM_EAST ccr
msmc2 -t 16 -i 20 -p $p -s -I 8-16,8-17,8-18,8-19,9-16,9-17,9-18,9-19,10-16,10-17,10-18,10-19,11-16,11-17,11-18,11-19 -o ${OUTDIR}/$dir/RM_EAST_cross_msmc.out $INPUTS

# WAC_EAST ccr
msmc2 -t 16 -i 20 -p $p -s -I 12-16,12-17,12-18,12-19,13-16,13-17,13-18,13-19,14-16,14-17,14-18,14-19,15-16,15-17,15-18,15-19 -o ${OUTDIR}/$dir/WAC_EAST_cross_msmc.out $INPUTS

# LAS_ORC ccr
msmc2 -t 16 -i 20 -p $p -s -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -o ${OUTDIR}/$dir/LAS_ORC_cross_msmc.out $INPUTS

# LAS_RM ccr
msmc2 -t 16 -i 20 -p $p -s -I 0-8,0-9,0-10,0-11,1-8,1-9,1-10,1-11,2-8,2-9,2-10,2-11,3-8,3-9,3-10,3-11 -o ${OUTDIR}/$dir/LAS_RM_cross_msmc.out $INPUTS

# LAS_WAC ccr
msmc2 -t 16 -i 20 -p $p -s -I 0-12,0-13,0-14,0-15,1-12,1-13,1-14,1-15,2-12,2-13,2-14,2-15,3-12,3-13,3-14,3-15 -o ${OUTDIR}/$dir/LAS_WAC_cross_msmc.out $INPUTS

# ORC_RM ccr
msmc2 -t 16 -i 20 -p $p -s -I 4-8,4-9,4-10,4-11,5-8,5-9,5-10,5-11,6-8,6-9,6-10,6-11,7-8,7-9,7-10,7-11 -o ${OUTDIR}/$dir/ORC_RM_cross_msmc.out $INPUTS

# ORC_WAC ccr
msmc2 -t 16 -i 20 -p $p -s -I 4-12,4-13,4-14,4-15,5-12,5-13,5-14,5-15,6-12,6-13,6-14,6-15,7-12,7-13,7-14,7-15 -o ${OUTDIR}/$dir/ORC_WAC_cross_msmc.out $INPUTS

# WAC_RM ccr
msmc2 -t 16 -i 20 -p $p -s -I 12-8,12-9,12-10,12-11,13-8,13-9,13-10,13-11,14-8,14-9,14-10,14-11,15-8,15-9,15-10,15-11 -o ${OUTDIR}/$dir/ORC_RM_cross_msmc.out $INPUTS


#### n = 1

## RUS

# SN_RM ccr
$msmc2 -t 16 -i 20 -p $p -s -I 24-20,24-21,25-20,25-21 -o ${OUTDIR}/$dir/SN_RUS_cross_msmc.out $INPUTS

# AK_RUS ccr
$msmc2 -t 16 -i 20 -p $p -s -I 22-20,22-21,23-20,23-21 -o ${OUTDIR}/$dir/AK_RUS_cross_msmc.out $INPUTS

# ECAN_RUS 
$msmc2 -t 16 -i 20 -p $p -s -I 16-20,16-21,17-20,17-21 -o ${OUTDIR}/$dir/ECAN_RUS_cross_msmc.out $INPUTS

# RM_RUS ccr
$msmc2 -t 16 -i 20 -p $p -s -I 8-20,8-20,9-20,9-21 -o ${OUTDIR}/$dir/RM_RUS_cross_msmc.out $INPUTS

# LAS_RUS ccr
$msmc2 -t 16 -i 20 -p $p -s -I 0-20,0-21,1-20,1-21 -o ${OUTDIR}/$dir/LAS_RUS_cross_msmc.out $INPUTS

# ORC_RUS ccr      
$msmc2 -t 16 -i 20 -p $p -s -I 4-20,4-21,5-20,5-21 -o ${OUTDIR}/$dir/ORC_RUS_cross_msmc.out $INPUTS

# WAS_RUS  
$msmc2 -t 16 -i 20 -p $p -s -I 12-20,12-21,13-20,13-21 -o ${OUTDIR}/$dir/WAS_RUS_cross_msmc.out $INPUTS

## AK

# SN_AK ccr
$msmc2 -t 16 -i 20 -p $p -s -I 24-22,24-23,25-22,25-23 -o ${OUTDIR}/$dir/SN_AK_cross_msmc.out $INPUTS

# ECAN_AK
$msmc2 -t 16 -i 20 -p $p -s -I 16-22,16-23,17-22,17-23 -o ${OUTDIR}/$dir/ECAN_AK_cross_msmc.out $INPUTS

# RM_AK ccr
$msmc2 -t 16 -i 20 -p $p -s -I 8-22,8-23,9-22,9-23 -o ${OUTDIR}/$dir/RM_AK_cross_msmc.out $INPUTS

# LAS_AK ccr
$msmc2 -t 16 -i 20 -p $p -s -I 0-22,0-23,1-22,1-23 -o ${OUTDIR}/$dir/LAS_AK_cross_msmc.out $INPUTS

# ORC_RM ccr      
$msmc2 -t 16 -i 20 -p $p -s -I 4-22,4-23,5-22,5-23 -o ${OUTDIR}/$dir/ORC_AK_cross_msmc.out $INPUTS

# WAS_AK  
$msmc2 -t 16 -i 20 -p $p -s -I 12-22,12-23,13-22,13-23 -o ${OUTDIR}/$dir/WAS_AK_cross_msmc.out $INPUTS

# LAS_AK ccr
$msmc2 -t 16 -i 20 -p $p -s -I 0-22,0-23,1-22,1-23 -o ${OUTDIR}/$dir/LAS_AK_cross_msmc.out $INPUTS

# ORC_RM ccr      
$msmc2 -t 16 -i 20 -p $p -s -I 4-22,4-23,5-22,5-23 -o ${OUTDIR}/$dir/ORC_AK_cross_msmc.out $INPUTS

# WAS_AK  
$msmc2 -t 16 -i 20 -p $p -s -I 12-22,12-23,13-22,13-23 -o ${OUTDIR}/$dir/WAS_AK_cross_msmc.out $INPUTS

##### combine into one file for R

~/bin/msmc-tools/combineCrossCoal.py LAS_EAST_msmc.out.final.txt LAS_msmc.out.final.txt EAST_msmc.out.final.txt > LAS_EAST_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py ORC_EAST_cross_msmc.out.final.txt ORC_msmc.out.final.txt EAST_msmc.out.final.txt > ORC_EAST_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py RM_EAST_cross_msmc.out.final.txt RM_msmc.out.final.txt EAST_msmc.out.final.txt > RM_EAST_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py WAC_EAST_cross_msmc.out.final.txt RM_msmc.out.final.txt EAST_msmc.out.final.txt > WAC_EAST_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py LAS_ORC_cross_msmc.out.final.txt LAS_msmc.out.final.txt ORC_msmc.out.final.txt > LAS_ORC_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py LAS_RM_cross_msmc.out.final.txt LAS_msmc.out.final.txt RM_msmc.out.final.txt > LAS_RM_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py LAS_WAC_cross_msmc.out.final.txt LAS_msmc.out.final.txt WAC_msmc.out.final.txt > LAS_WAC_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py ORC_RM_cross_msmc.out.final.txt ORC_msmc.out.final.txt RM_msmc.out.final.txt > ORC_RM_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py ORC_WAC_cross_msmc.out.final.txt ORC_msmc.out.final.txt WAC_msmc.out.final.txt > ORC_WAC_msmc.combined.final.txt

~/bin/msmc-tools/combineCrossCoal.py WAC_RM_cross_msmc.out.final.txt WAC_msmc.out.final.txt RM_msmc.out.final.txt > WAC_RM_msmc.combined.final.txt

# this part was run on ceres, had to change paths and module/conda info

export PATH="/home/catherine.quinn/miniforge3/bin:$PATH"
source /home/catherine.quinn/miniforge3/etc/profile.d/conda.sh 
source /home/catherine.quinn/miniforge3/etc/profile.d/mamba.sh
source activate /home/catherine.quinn/miniforge3/envs/demographic
mamba activate demographic
msmc2="msmc2_Linux"

combine="/home/catherine.quinn/bin/msmc-tools/combineCrossCoal.py"

cd /project/conservationgenomicsquinn/redfoxes/msmc2/CCR

$combine SN_RUS_cross_msmc.out.final.txt SN_msmc.out.final.txt RUS_msmc.out.final.txt > SN_RUS_msmc.combined.final.txt
$combine AK_RUS_cross_msmc.out.final.txt AK_msmc.out.final.txt RUS_msmc.out.final.txt > AK_RUS_msmc.combined.final.txt
$combine ECAN_RUS_cross_msmc.out.final.txt ECAN_S14-0358_msmc.out.final.txt RUS_msmc.out.final.txt > ECAN_RUS_msmc.combined.final.txt
$combine RM_RUS_cross_msmc.out.final.txt RM_msmc.out.final.txt RUS_msmc.out.final.txt > RM_RUS_msmc.combined.final.txt

$combine LAS_RUS_cross_msmc.out.final.txt LAS_msmc.out.final.txt RUS_msmc.out.final.txt > LAS_RUS_msmc.combined.final.txt
$combine ORC_RUS_cross_msmc.out.final.txt ORC_msmc.out.final.txt RUS_msmc.out.final.txt > ORC_RUS_msmc.combined.final.txt
$combine WAS_RUS_cross_msmc.out.final.txt WAC_msmc.out.final.txt RUS_msmc.out.final.txt > WAS_RUS_msmc.combined.final.txt


$combine SN_AK_cross_msmc.out.final.txt SN_msmc.out.final.txt AK_msmc.out.final.txt > SN_AK_msmc.combined.final.txt
$combine ECAN_AK_cross_msmc.out.final.txt ECAN_S14-0358_msmc.out.final.txt AK_msmc.out.final.txt > ECAN_AK_msmc.combined.final.txt
$combine RM_AK_cross_msmc.out.final.txt RM_msmc.out.final.txt AK_msmc.out.final.txt > RM_AK_msmc.combined.final.txt

$combine LAS_AK_cross_msmc.out.final.txt LAS_msmc.out.final.txt AK_msmc.out.final.txt > LAS_AK_msmc.combined.final.txt
$combine ORC_AK_cross_msmc.out.final.txt ORC_msmc.out.final.txt AK_msmc.out.final.txt > ORC_AK_msmc.combined.final.txt
$combine WAS_AK_cross_msmc.out.final.txt WAC_msmc.out.final.txt AK_msmc.out.final.txt > WAS_AK_msmc.combined.final.txt
