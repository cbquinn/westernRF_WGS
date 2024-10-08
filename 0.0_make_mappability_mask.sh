#!/bin/bash -l
#SBATCH --job-name=gem
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --time 10-00:00:00
#SBATCH --mem=64GB
#SBATCH -p bmm
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/make_mappability_mask_gem.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/make_mappability_mask_gem.err

# Make mappability mask with [genmap](https://github.com/cpockrandt/genmap)

# make things fail on errors
set -o nounset
set -o errexit

dir="/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1"
ref="${dir}/GCF_018345385.1_ASM1834538v1_genomic.fna"
out="$dir/genmap/artic"
scaff=GCF_018345385.1_ASM1834538v1_genomic
autosomes=/group/ctbrowngrp2/cbquinn/fox4/ref/GCF_018345385.1/autosomes.bed

start=`date +%s`
echo $HOSTNAME

module load bedtools2
module load repeatmasker
conda activate genmap

# echo "Indexing fasta..."
genmap index -S 45 -F $ref -I $out

#echo "computing mapability..."
genmap map -K 100 -E 2 -I ${out} -O ${out} -bg

#echo "Finding repeat regions..."

RepeatMasker -noint $ref -dir ${out}

#echo "Combining results into single bed file..."

# convert RepeatMasker output into bed file
#tail -n +4 ${out}/${scaff}.fa.out | sed -E 's/.+(SUPER__.+)/\1/' | #sed -E 's/\s+/\t/g' | awk '{print $1"\t"$2-1"\t"$3-1}' > ${out}/${scaff}_repeat_mask.bed

# convert genmap mappability output into bed file
awk '{if ($4 != 1) print $0}' ${out}/${scaff}.genmap.bedgraph > ${out}/${scaff}_E2_low_mapability_mask.bed

# convert genmap mappability mask to good scores to keep
bedtools intersect \
	-a ${out}/${scaff}_E2_low_mapability_mask.bed \
	-b $autosomes \
	-wa | \
bedtools complement \
	-i stdin \
	-g <(cut -f1,3 $autosomes) > ${out}/${scaff}_E2_high_mapability_mask.bed

end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
