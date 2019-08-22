#!/bin/bash
#$ -cwd

declare -a traits=("cad_cases" "stroke_cases" "copd_cases"
                "dementia_cases" "lungcancer_cases" "colorectal_cancer_cases"
                "renal_failure_cases" "breast_cancer_cases"
                "t2d_cases_prob" "t2d_cases_probposs" "t1d_cases_prob"
                "t1d_cases_probposs" "any_infertility_cases"
                "nafld_cases" "cld_cases" "smoker_cases" "aki_cases" "ckd_cases" "haem_stroke_cases" "isch_stroke_cases")
declare -a eth_groups=("all.white")
declare -a sex_groups=("comb" "men" "women")

#Collinearity between age_assessment and age_sqaured 
#and for some traits perfect pseudoseparation of cases if include assessment centre, 
#thus only adjust for the things below
for trait in "${traits[@]}"
do
        for eth_group in "${eth_groups[@]}"
        do
		for sex_group in "${sex_groups[@]}"
		do
                $plink2 --bfile ../pulit.180916 \
                --keep ../sex.het.enrichment/${eth_group}.${sex_group}.samples.to.keep \
                --threads 7 \
                --memory 105000 \
                --glm no-x-sex \
                --pheno ../sex.het.enrichment/plink.phenotype.file.from.ukbb.cases.for.logistic.180921.txt \
                --pheno-name ${trait} \
                --1 \
                --covar ../sex.het.enrichment/plink.phenotype.file.from.ukbb.cases.for.logistic.180921.txt \
                --covar-name age_assessment PC1-PC10 dummy_array dummy_sex \
                --out ../sex.het.enrichment/${eth_group}.${sex_group}.pulit
		done
        done
done




