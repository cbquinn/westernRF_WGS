
# get q values from run with highest likelihood
# get CV values from all runs
# (see https://github.com/dloesch/leverage-ancestry)

INDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/final
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/admixture/final/all
data="$INDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked_LDpruned.fam"
# or

INDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/final
OUTDIR=/group/ctbrowngrp2/cbquinn/fox4/1_vcfs/variant/n34/admixture/final/reduced
data="$INDIR/bcf_variant_filterA_ac1.dp3.mis20_n34_masked_LDpruned_reduced.fam"

cd $OUTDIR

for K in {1..15}
do
cd K${K}

echo "CV_error" > K${K}.CV_error.txt
echo "loglikelihood" > K${K}.loglikelihood.txt

for run in {1..10}
do
#CV error
cat ./run${run}/run${run}.K${K}.o* | grep CV >> K${K}.CV_error.txt

#loglikelihood
cat ./run${run}/run${run}.K${K}.o* | grep Log | tail -1 >> K${K}.loglikelihood.txt
done

#best loglikelihood
best=$(cat K${K}.loglikelihood.txt | grep Log | cut -d" " -f2 | sort -n | tail -1 | sed 's/-//g')

#summary file
paste K${K}.random.txt K${K}.CV_error.txt K${K}.loglikelihood.txt > K${K}.summary.txt

#best run
cat K${K}.summary.txt | grep "$best" > K${K}.best_run.txt

#create outputfile from best run with column ids
#column id files were made from IID column of PED file

cat $data | cut -f1 -d" " > K${K}.subjects.txt
 
best=$(cat K${K}.best_run.txt | sed "1q;d" | cut -f1 -d" ")
results=$(ls ./$best | grep Q)
paste K${K}.subjects.txt ./$best/$results > K${K}.results.txt

cd ..
done

# get CV errors for best runs
rm bestKruns_statssummary.txt
for K in {1..15}
do
cat K${K}/K${K}.best_run.txt >> bestKruns_statssummary.txt
done
