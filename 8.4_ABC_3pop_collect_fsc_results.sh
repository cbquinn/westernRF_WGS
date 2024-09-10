#!/bin/bash -l
#SBATCH --job-name=jsfs
#SBATCH --time 04:00:00
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH -p bml
#SBATCH -A ctbrowngrp
#SBATCH -o /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/fastsimcoal_jsfs.out
#SBATCH -e /group/ctbrowngrp2/cbquinn/fox4/slurmlogs/fastsimcoal_jsfs.err

# convert from 2d to 1d array and remove the monomorphic reference bin
# but for meta, we need to do this 3 times


cd /group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240719meta/out


s=${1}

# this is for each array run
for run in {1..50}
    do

        # there are 3 jsfs generated per iteration
        for popcombo in pop1_0 pop2_0 pop2_1
        do

            rm ${s}_${run}/jsfs_${s}_${run}_${popcombo}.txt
            echo ${s}_${run}_${popcombo}
                
                # there are 10000 iterations for each
                for i in {1..1000}
                do
                    FILE=${s}_${run}/${s}_${run}_def${i}_jointDAF${popcombo}.obs
                    # Check if the file exists
                        if [ -e "$FILE" ]; then
                                    sed '1,2d' ${FILE} | cut -f2- | tr '\n' '\t' | cut -f2- >> ${s}_${run}/jsfs_${s}_${run}_${popcombo}.txt
                            else
                                    echo "File not found: $FILE"
                                    break  # Exit the inner loop if the file is not present
                    fi
                done
            done
        done
        
