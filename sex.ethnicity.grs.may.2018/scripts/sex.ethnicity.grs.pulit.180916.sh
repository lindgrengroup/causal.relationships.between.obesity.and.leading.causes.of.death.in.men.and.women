#!/bin/bash
#$ -cwd

##############################################################
######## Make sex- and ethnicity obesity GRS on UKBB data ########
############### Jenny Censin, 2018-05-13 #####################

#Settings
declare -a traits=("bmi" "whr" "whradjbmi")
declare -a snp_groups=("comb.pulit.sig" "men.pulit"
		"men.pulit.sig" "men.pulit.phet"
		"women.pulit" "women.pulit.sig"
		"women.pulit.phet"
		"comb.pulit.winner" "men.pulit.winner" "women.pulit.winner"
		"comb.pulit.winner_unweighted")

plink=/apps/well/plink/1.90b3/plink
export plink

for trait in "${traits[@]}"
do
	for snp_group in "${snp_groups[@]}"
	do
	$plink --bfile ../pulit.180916 \
	--score ../${trait}.eur.${snp_group}.180513.txt 1 4 8 header sum \
        --threads 1 \
	--memory 15000 \
        --out ../bmi.whr.profiles/grs.${trait}.eur.${snp_group}.180513
	done
done


