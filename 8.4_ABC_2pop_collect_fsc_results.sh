#!/bin/bash -l
#SBATCH --job-name=jsfs
#SBATCH --time 04:00:00
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH -p bml
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/fastsimcoal_jsfs.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/fastsimcoal_jsfs.err

# takes sfs created from fsc and...
# converts them from 2d to 1d array (all entries on a single line)
# removes the monomorphic reference bin
# collects them into a single file

cd /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240625/out

s=${1}

for run in {1..10}
        do
        rm ${s}_${run}/jsfs_${s}_${run}.txt
        echo ${s}_${run}
        for i in {1..10000}
        do
        FILE=${s}_${run}/${s}_${run}_def${i}_jointDAFpop1_0.obs
        # Check if the file exists
            if [ -e "$FILE" ]; then
                        sed '1,2d' ${FILE} | cut -f2- | tr '\n' '\t' | cut -f2- >> ${s}_${run}/jsfs_${s}_${run}.txt
                else
                        echo "File not found: $FILE"
                        break  # Exit the inner loop if the file is not present
            fi
        done
        done
