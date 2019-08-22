#!/bin/bash
#$ -cwd

#Settings
declare -a traits=("bmi" "whr" "whradjbmi")
declare -a sex_groups=("comb" "men" "women")
declare -a snp_groups=("normal" "half" "randomnoise" "systematicnoise" "unweighted")

plink=/apps/well/plink/1.90b3/plink
export plink

for trait in "${traits[@]}"
do
	for sex_group in "${sex_groups[@]}"
	do 
		for snp_group in "${snp_groups[@]}"
	        do

	        $plink --bfile ../pulit.180916 \
        	--score ../sensitivity.risk.weights/${trait}.eur.${sex_group}.internal.${snp_group}.180513.txt 1 2 8 header sum \
	        --threads 1 \
        	--memory 15000 \
	        --out ../bmi.whr.profiles/grs.${trait}.eur.${sex_group}.internal.${snp_group}.180513
        	done
	done
done

declare -a snp_groups=("giukbb")

for trait in "${traits[@]}"
do
        for sex_group in "${sex_groups[@]}"
        do
                for snp_group in "${snp_groups[@]}"
                do

                $plink --bfile ../pulit.180916 \
                --score ../sensitivity.risk.weights/${trait}.eur.${sex_group}.${snp_group}.sig.180513.txt 1 4 7 header sum \
                --threads 1 \
                --memory 15000 \
                --out ../bmi.whr.profiles/grs.${trait}.eur.${sex_group}.${snp_group}.sig.180513
                done
        done
done

declare -a snp_groups=("0.01.fdr" "0.05.fdr" "0.1.fdr")
declare -a sex_groups=("men" "women")
for trait in "${traits[@]}"
do
        for sex_group in "${sex_groups[@]}"
        do
                for snp_group in "${snp_groups[@]}"
                do

                $plink --bfile ../pulit.180916 \
                --score ../sensitivity.risk.weights/${trait}.eur.${sex_group}.${snp_group}.180513.txt 1 4 8 header sum \
                --threads 1 \
                --memory 15000 \
                --out ../bmi.whr.profiles/grs.${trait}.eur.${sex_group}.${snp_group}.180513
                done
        done
done

