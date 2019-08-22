#!/bin/env Rscript
#$ -cwd

#####################
##################### SCRIPT FOR MAKING SMALLER rs_aliases files

#Script for making the rsalias file and README
#File with all the rsAliases and how they've been updated in the NCBI database.
#File downloaded from: 
#ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/data/organism_data/
#Date downloaded: 2018-03-15
#Date last updated: 2017-03-13
#Original file name: RsMergeArch.bcp.gz

#Column names according to the file downloaded from here: 
#ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/data/organism_schema/
#Date downloaded: 2018-03-15
#Original file name: human_9606_table
#Original column names according to the file, the actual file don't have any column 
#names though:  
#[rsHigh] [bigint] NOT NULL ,
#[rsLow] [bigint] NOT NULL ,
#[build_id] [int] NULL ,
#[orien] [tinyint] NOT NULL ,
#[create_time] [datetime] NOT NULL ,
#[last_updated_time] [datetime] NOT NULL ,
#[rsCurrent] [int] NULL ,
#[orien2Current] [tinyint] NULL ,
#[comment] [varchar](255) NULL

################## OLD ONE!

#Read in table and give correct column names:
df <- read.table("../RsMergeArch.bcp.gz", stringsAsFactors = F, header = F, sep = "\t")
colnames(df) <- c("rsHigh", "rsLow", "buildID", "orien", "create_time", "updated_time", "rsCurrent", "orien2current", "comment")

#Subset the relevant columns and update so that the SNPs have "rs" in front of number:
df <- df[, c(1,7)]
df$rsHigh <- paste("rs", df$rsHigh, sep = "")
df$rsCurrent <- paste("rs", df$rsCurrent, sep = "")
#Write the table to be used:
write.table(df, "../rsaliases.txt", quote = F, row.names = F)

################### NEW ONE - see README file for where file is from: 
#The file RsMergeArch.180920.txt.gz was downloaded from
#ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/
#File name: RsMergeArch.bcp.gz
#Date downloaded: 2018-09-20
#Date last uploaded according to website: 2018-02-07
#Checked the file

#Read in the file
df <- read.table(gzfile("../RsMergeArch.180920.txt.gz"), stringsAsFactors = F, header = F, sep = "\t")

#50 SNPs (the first ones) are NA in rsCurrent. The rsHigh of those aren't online on dbSNPs website, 
#however the rsLow is. So, update to that. 
df[is.na(df$V7) & df$V3 == 151, "V7"] <- df[is.na(df$V7) & df$V3 == 151, "V2"]

df <- df[, c(1,7)]
colnames(df) <- c("rsHigh", "rsCurrent")
df$rsHigh <- paste("rs", df$rsHigh, sep = "")
df$rsCurrent <- paste("rs", df$rsCurrent, sep = "")

#Write the table to be used:
write.table(df, "../rsaliases.181010.txt", quote = F, row.names = F)

