#!/bin/bash
#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 7
#$ -N mr.cont
#$ -o ../oldjobs
#$ -j y

echo "########################################################"
echo "GRS submit mr cont"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

module load R/3.4.3
#Pathway to the script that does the MRs
R --vanilla --file=./mr.adiposity.traits.to.continuous.traits.ipd.R

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

