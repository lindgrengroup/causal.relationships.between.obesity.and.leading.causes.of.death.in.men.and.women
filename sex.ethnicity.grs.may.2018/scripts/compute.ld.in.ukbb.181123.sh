#!/bin/bash
#$-cwd

#mkdir temp.ld.between.snps.181123

#Lists of SNPs for each trait: 
declare -a traits=("bmi" "whr" "whradjbmi")
declare -a sex_groups=("men" "women")

for trait in "${traits[@]}"
do
	tail -n +2 ../index.files.with.trialleles/${trait}.eur.comb.pulit.sig.180513.txt | 
		cut -d " " -f1 > ../temp.ld.between.snps.181123/${trait}.comb.pulit.snps

	tail -n +2 ../index.files.with.trialleles/${trait}.eur.comb.pulit.winner.180513.txt |
                cut -d " " -f1 > ../temp.ld.between.snps.181123/${trait}.comb.pulit_winner.snps

	for sex_group in "${sex_groups[@]}"
	do
		#Pulit
		tail -n +2 ../index.files.with.trialleles/${trait}.eur.${sex_group}.pulit.180513.txt | \
		cut -d " " -f1 > ../temp.ld.between.snps.181123/${trait}.${sex_group}.pulit.snps
	done
done

declare -a sex_groups=("comb" "men" "women")
declare -a datasets=("pulit" "pulit_winner")

for trait in "${traits[@]}"
do
	for sex_group in "${sex_groups[@]}"
	do
		for dataset in "${datasets[@]}"
		do
			if [ ${dataset} == "pulit" ] || ([ ${dataset} == "pulit_winner" ] && [ ${sex_group} == "comb" ])
			then
				$plink --bfile ../pulit.180916 \
				--keep /well/lindgren/jc/ukbb/ukbb.completely.unrelated.european.samples.190215.txt \
				--extract ../temp.ld.between.snps.181123/${trait}.${sex_group}.${dataset}.snps \
				--r2 inter-chr \
				--ld-window-r2 0 \
				--out ../temp.ld.between.snps.181123/${trait}.${sex_group}.${dataset}.snps.ld.181123
			fi
		done
	done
done


