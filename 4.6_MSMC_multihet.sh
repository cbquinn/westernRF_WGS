
#!/bin/bash -l
#SBATCH --job-name=multihet
#SBATCH --nodes 1
#SBATCH --array=1-24
#SBATCH --ntasks 1
#SBATCH --time 00:30:00
#SBATCH --mem=1GB
#SBATCH -p bmh
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_multihet_phased_%a.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/msmc2_multihet_phased_%a.err

CHR=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/autosomes.bed | cut -f1)
INDIR=/group/ctbrowngrp2/cbquinn/fox4/4_SMC/MSMC2/inputs
BED=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed

cd $INDIR

# masks are regions that can be kept because of mappability or sample-specific depth
~/bin/msmc-tools/generate_multihetsep.py \
--mask=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/genmap/artic/GCF_018345385.1_ASM1834538v1_genomic_E2_high_mapability_mask.bed \
--mask=LAS_F01_${CHR}_DEPTHmask.bed.gz \
--mask=LAS_F02_${CHR}_DEPTHmask.bed.gz \
--mask=ORC_S18-2071_${CHR}_DEPTHmask.bed.gz \
--mask=ORC_S17-2544_${CHR}_DEPTHmask.bed.gz \
--mask=RM_S13-2269_${CHR}_DEPTHmask.bed.gz \
--mask=RM_S13-2260_${CHR}_DEPTHmask.bed.gz \
--mask=WAC_S10-0511_${CHR}_DEPTHmask.bed.gz \
--mask=WAC_S11-0716_${CHR}_DEPTHmask.bed.gz \
--mask=ECAN_S14-0358_${CHR}_DEPTHmask.bed.gz \
--mask=VT_F12-232_${CHR}_DEPTHmask.bed.gz \
--mask=RUS_S12-0237_${CHR}_DEPTHmask.bed.gz \
--mask=AK_S12-1159_${CHR}_DEPTHmask.bed.gz \
--mask=SN_S11-0008_${CHR}_DEPTHmask.bed.gz \
phased/LAS_F01_${CHR}_phased_final.vcf.gz \
phased/LAS_F02_${CHR}_phased_final.vcf.gz \
phased/ORC_S18-2071_${CHR}_phased_final.vcf.gz \
phased/ORC_S17-2544_${CHR}_phased_final.vcf.gz \
phased/RM_S13-2269_${CHR}_phased_final.vcf.gz \
phased/RM_S13-2260_${CHR}_phased_final.vcf.gz \
phased/WAC_S10-0511_${CHR}_phased_final.vcf.gz \
phased/WAC_S11-0716_${CHR}_phased_final.vcf.gz \
phased/ECAN_S14-0358_${CHR}_phased_final.vcf.gz \
phased/VT_F12-232_${CHR}_phased_final.vcf.gz \
phased/RUS_S12-0237_${CHR}_phased_final.vcf.gz \
phased/AK_S12-1159_${CHR}_phased_final.vcf.gz \
phased/SN_S11-0008_${CHR}_phased_final.vcf.gz \
> phased/all13_${CHR}_msmc.txt
