#!/bin/bash
#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N grs.180911
#$ -o ../
#$ -j y
#$ -t 1-22

echo "########################################################"
echo "Make UKBB files to perform GRSs on"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Working directory: "`pwd`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

##############################################################
### Make bed files to perform the GRSs on in the UKBB    #####
############### Jenny Censin, 2018-09-11 #####################

plink2=/apps/well/plink/2.00a-20170724/plink2
export plink2

#Pathways to the UKBB imputed bgen files, the corresponding sample file, and previously made files
$plink2 --bgen /well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_imp_chr${SGE_TASK_ID}_v3.bgen \
--sample /well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--keep ../temporary.bed.pulit.180911/samples.to.keep.txt \
--threads 1 \
--extract ../temporary.bed.pulit.180911/snps.to.extract.txt \
--out ../temporary.bed.pulit.180911/grs.chr${SGE_TASK_ID} \
--memory 15000 \
--make-bed

if [ -f ../temporary.bed.pulit.180911/grs.chr${SGE_TASK_ID}.bed ]
then
echo "../temporary.bed.pulit.180911/grs.chr${SGE_TASK_ID}">>../temporary.bed.pulit.180911/files_to_merge.txt
fi

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0



