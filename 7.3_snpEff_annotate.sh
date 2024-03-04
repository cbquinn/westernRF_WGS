#!/bin/bash -l
#SBATCH --job-name=snpeff
#SBATCH --time 24:00:00
#SBATCH --mem=16GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/run_snpeff.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/run_snpeff.err

conda activate snpeff
module load tabix

VCF=/group/ctbrowngrp2/cbquinn/fox4/5_load/polarize/filterA_all_wGERP_AArotated.vcf.gz
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/5_load/snpeff
OUTFILE=filterA_all_wGERP_AArotated_snpeff_all
database=Vulpes_lagopus

cd $OUTDIR

#snpEff -Xmx14g -v \
#	-nodownload \
#	-csvStats ${OUTFILE}.csvStats.csv \
#	${database} ${VCF} | bgzip -c > ${OUTFILE}.vcf.gz

snpEff -Xmx14g -v \
	-nodownload \
	${database} ${VCF} | bgzip -c > ${OUTFILE}.vcf.gz

bcftools index ${OUT}.vcf.gz
```

