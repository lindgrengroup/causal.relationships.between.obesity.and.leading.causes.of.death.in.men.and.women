#!/bin/bash
#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -N submit.sample.qc
#$ -o /well/lindgren/jc/ukbb/
#$ -j y

echo "########################################################"
echo "Submit sample QC job"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

module load R/3.4.3
#Pathway to the script that does sample QC
R --vanilla --file=./sample.qc.ukbb.all.ethnicities.180507.R  

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

