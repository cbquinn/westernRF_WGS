#!/bin/bash -l
#SBATCH --job-name=bcfroh
#SBATCH --time 1-00:00:00
#SBATCH --mem=6GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/bcftools_ROH.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/bcftools_ROH.err

# run bcftools/roh with default allele frequency and mean recombinatino 

cd /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/3_ROH/bcftools

module load bcftools/1.14

# this creates output files for filters A,B,C,D (used A in final analysis)

for i in $(ls *_masked.vcf.gz)
do
start=`date +%s`
bcftools roh --G 30 --rec-rate 1e-6 --AF-dflt 0.4 -o $OUTDIR/${i%.vcf.gz}.bcf.ROH $i
grep ^RG $OUTDIR/${i%.vcf.gz}.bcf.ROH > $OUTDIR/${i%.vcf.gz}.bcf.ROH_RG.txt
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME for $i: $runtime"
done
