#!/bin/bash -l
#SBATCH --job-name=plinkroh
#SBATCH --time 1-00:00:00
#SBATCH --mem=3GB
[[3.3_IBS_ROH]]#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/plink_ROH.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/plink_ROH.err

# identify ROH with Plink rule-based algorithm for comparison with bcftools/ROH
# also ran for all 5 filters, and testing 5 distinct combos of rules

cd /group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/3_ROH/plink
INDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34

module load deprecated/plink/1.90

cd $INDIR
for i in $(ls *_masked.vcf.gz)
do
VCF=$i
OUT=${i%.vcf.gz}

# create plink files
plink --vcf $VCF --keep-allele-order --allow-extra-chr --chr-set 24 --recode --out $INDIR/$OUT

plink --file $INDIR/$OUT \
    --allow-extra-chr \
    --chr-set 24 \
    --homozyg-window-snp 50 \
    --homozyg-window-het 3 \
    --homozyg-window-missing 10 \
    --homozyg-window-threshold 0.05 \
    --homozyg-snp 50 \
    --homozyg-kb 100 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --out $OUTDIR/${OUT}_plink_winsnp50.winhet3.winmis10.winthresh10.hom50.homkb100.homdens50.homgap1000

plink --file $INDIR/$OUT \
    --allow-extra-chr \
    --chr-set 24 \
    --homozyg-window-snp 50 \
    --homozyg-window-het 3 \
    --homozyg-window-missing 10 \
    --homozyg-window-threshold 0.05 \
    --homozyg-snp 25 \
    --homozyg-kb 100 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --out $OUTDIR/${OUT}_plink_winsnp25.winhet3.winmis10.winthresh10.hom50.homkb100.homdens50.homgap1000

plink --file $INDIR/$OUT \
    --allow-extra-chr \
    --chr-set 24 \
    --homozyg-window-snp 100 \
    --homozyg-window-het 1 \
    --homozyg-window-missing 20 \
    --homozyg-window-threshold 0.05 \
    --homozyg-snp 25 \
    --homozyg-kb 100 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --homozyg-het 750 \
    --out $OUTDIR/${OUT}_plink_winsnp100.winhet1.winmis20.winthresh05.hom25.homkb100.homdens50.homgap1000.homhet750
    
plink --file $INDIR/$OUT \
    --allow-extra-chr \
    --chr-set 24 \
    --homozyg-window-snp 250 \
    --homozyg-window-het 3 \
    --homozyg-window-missing 20 \
    --homozyg-window-threshold 0.05 \
    --homozyg-snp 25 \
    --homozyg-kb 100 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --homozyg-het 750 \
    --out $OUTDIR/${OUT}_plink_winsnp250.winhet3.winmis20.winthresh05.hom25.homkb100.homdens50.homgap1000.homhet750

plink --file $INDIR/$OUT \
    --allow-extra-chr \
    --chr-set 24 \
    --homozyg-window-snp 500 \
    --homozyg-window-het 5 \
    --homozyg-window-missing 20 \
    --homozyg-window-threshold 0.05 \
    --homozyg-snp 25 \
    --homozyg-kb 100 \
    --homozyg-density 50 \
    --homozyg-gap 1000 \
    --homozyg-het 750 \
    --out $OUTDIR/${OUT}_plink_winsnp500.winhet5.winmis20.winthresh05.hom25.homkb100.homdens50.homgap1000.homhet750
    
done
