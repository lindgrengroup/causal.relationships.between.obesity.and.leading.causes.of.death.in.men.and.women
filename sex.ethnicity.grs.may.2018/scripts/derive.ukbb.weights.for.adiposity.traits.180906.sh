#!/bin/bash
#$ -cwd

plink=/apps/well/plink/1.90b3/plink
export plink

#Create list of individuals to keep for all analyses
awk 'BEGIN {FS = "[ \t]+"} \
{ if (NR != 1) { print $1 " " $1 }}' ../ukbb.residuals.for.adiposity.traits.180906.txt \
> ../temporary.ukbb.weights/comb

awk 'BEGIN {FS = "[ \t]+"} \
{ if ($3 == "M") print $1 " " $1 }' ../ukbb.residuals.for.adiposity.traits.180906.txt \
> ../temporary.ukbb.weights/men

awk 'BEGIN {FS = "[ \t]+"} \
{ if ($3 == "F") print $1 " " $1}' ../ukbb.residuals.for.adiposity.traits.180906.txt \
> ../temporary.ukbb.weights/women

declare -a traits=("bmi" "whr" "res_whr_inv")
declare -a sex_groups=("comb" "women" "men")

for trait in "${traits[@]}"
do
	for sex_group in "${sex_groups[@]}"
	do
		$plink --bfile ../pulit.180916 \
		--keep ../temporary.ukbb.weights/${sex_group} \
		--threads 7 \
		--memory 105000 \
		--linear \
		--pheno ../ukbb.residuals.for.adiposity.traits.180906.txt \
		--pheno-name ${trait}_${sex_group} \
		--out ../temporary.ukbb.weights/${trait}_${sex_group}
	done
done

