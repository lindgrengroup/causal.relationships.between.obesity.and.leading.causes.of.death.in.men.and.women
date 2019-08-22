#!/bin/bash
#$ -cwd

#Make output directories
#mkdir ../temporary.ukbb.bed.180507
#mkdir ../genotype.qc.snplists.180507

#Make list of the individuals that pass QC, non-Europeans included
#the pathway is to the file previously made with all samples passing QC
tail -n +2 ../ukbb.samples.passing.qc.relevant.pheno.all.ancestries.180504.txt | \
cut -d ' ' -f1 > ../temporary.ukbb.bed.180507/samples.to.keep.txt

#Make a file with unrelated European individuals to compute
#QC in for projects with just Europeans
awk -F" " '{if ($57=="brit" || $57=="white" || $57 == "recode.white") \
print $1}' ../ukbb.samples.passing.qc.relevant.pheno.all.ancestries.180504.txt \
> ../temporary.ukbb.bed.180507/euro.samples.to.keep.txt

make_bed_files() {
	#Creates a temporary BED file for each chromosome. The bgen pathway is to the 
	#UKBB imputation files, the sample file the corresponding sample file. 
	$plink2 --bgen /well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_imp_chr${SGE_TASK_ID}_v3.bgen \
	--sample /well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
	--keep ../temporary.ukbb.bed.180507/euro.samples.to.keep.txt \
	--threads 1 \
	--memory 15000 \
	--hard-call-threshold 0.1 \
	--make-bed \
	--out ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}

	#Makes the SNPs into the rs123_A_C format
	awk 'BEGIN {FS = "[[:space:]]+"} \
	{a[$5]; a[$6]; asorti(a,b)
	$2=$2"_"b[1]"_"b[2]; delete a}1' \
	../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}.bim > \
	../temporary.ukbb.bed.180507/first.${SGE_TASK_ID}

	rm ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}.bim
	mv ../temporary.ukbb.bed.180507/first.${SGE_TASK_ID} ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}.bim

	#Get list of SNPs that pass geno
	$plink --bfile ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID} \
	--geno 0.05 \
	--threads 1 \
	--memory 15000 \
	--out ../genotype.qc.snplists.180507/geno.180507.chr${SGE_TASK_ID} \
	--write-snplist

	#Get lists of SNPs that pass HWE
	$plink --bfile ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID} \
	--hwe 0.000001 \
	--threads 1 \
	--memory 15000 \
	--out ../genotype.qc.snplists.180507/hwe.180507.chr${SGE_TASK_ID} \	
	--write-snplist

	#Get lists of SNPs that pass MAF
	$plink --bfile ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID} \
	--maf 0.0001 \
	--threads 1 \
	--memory 15000 \
	--out ../genotype.qc.snplists.180507/maf.180507.chr${SGE_TASK_ID} \
	--write-snplist

	#Get the SNPs that have an INFO score > 0.3 and make into rs123_A_C format, 
	#with the pathway being to the UKBB provided genotype QC files
	awk 'BEGIN {FS = "[[:space:]]+"} \	
	{a[$4]; a[$5]; asorti(a,b)
	$2=$2"_"b[1]"_"b[2]; delete a}
	{ if ($8>0.3 && $8!="NA") print $2 }' \
	/well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_mfi_chr${SGE_TASK_ID}_v3.txt \
	> ../genotype.qc.snplists.180507/info.180507.chr${SGE_TASK_ID}

	rm ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}.bed
}

for SGE_TASK_ID in {1..22}
do
        make_bed_files "SGE_TASK_ID" &
done
wait

for SGE_TASK_ID in {1..22}
do
	#Identify duplicated SNP positions by position
	cat ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}.bim | cut -d" " -f4 |
	sort | uniq -c | awk '$1>1{print $2}' > ../temporary.ukbb.bed.180507/pos_non_biallelic_snps_chr${SGE_TASK_ID}.181025.txt

        #Find duplicated SNPs from the SNPs that pass all QC by position, and make into unique IDs
        awk 'NR==FNR{a[$0];next}($4 in a)' ../temporary.ukbb.bed.180507/pos_non_biallelic_snps_chr${SGE_TASK_ID}.181025.txt \
        ../temporary.ukbb.bed.180507/raw.ukbb.sample.qc.${SGE_TASK_ID}.bim > ../temporary.ukbb.bed.180507/non_biallelic_snps_chr${SGE_TASK_ID}.txt

        #Write list of non-biallelic SNPs
        cat ../temporary.ukbb.bed.180507/non_biallelic_snps_chr${SGE_TASK_ID}.txt | cut -d" " -f2 >> ../genotype.qc.snplists.180507/non.biallelic.snps.txt
done


