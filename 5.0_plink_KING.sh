#!/bin/bash -l
#SBATCH --job-name=king
#SBATCH --time 08:00:00
#SBATCH --mem=4GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/plink_KING.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/plink_KING.err


OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/KING
INDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34

module load deprecated/plink2/2.3-alpha-final

cd $INDIR
for i in $(ls *_masked.vcf.gz)
do
VCF=$i
OUT=${i%.vcf.gz}

# create plink input files
plink2 --vcf $VCF \
    --allow-extra-chr \
    --chr-set 24 \
    --make-bed \
    --out ${OUT}
    
# run KING
plink2 --bfile $INDIR/$OUT \
    --allow-extra-chr \
    --chr-set 24 \
	--make-king 'square' --out KING/${OUT}
done

