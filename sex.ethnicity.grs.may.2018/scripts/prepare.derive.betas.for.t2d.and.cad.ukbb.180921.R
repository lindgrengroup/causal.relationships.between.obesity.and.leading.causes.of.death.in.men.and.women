#!/bin/env Rscript
#$ -cwd

df <- read.table("../ukbb_cases_for_logistic_regression_180620.txt", stringsAsFactors = F, 
	header = T)

#Make dummy variables for assessment centre
new_rows <- sapply(unique(unlist(df$assessment_centre)), function(x) as.numeric(grepl(x, df$assessment_centre)))
colnames(new_rows) <- paste("assessment_centre_", unique(unlist(df$assessment_centre)), sep = "")

df <- cbind(df, new_rows)

df$dummy_sex <- ifelse(df$Submitted_Gender == "M", 0, ifelse(df$Submitted_Gender == "F", 1, "MISSING"))
df$dummy_array <- ifelse(df$genotyping_array == "UKBB", 0, ifelse(df$genotyping_array == "UKBL", 1, "MISSING"))

#Change the ones to exclude to -9 so works with plink
case_columns <- colnames(df)[grep("_cases", colnames(df))]
df[, case_columns] <- apply(df[, case_columns], 2, function(x) ifelse(x == 2, -9, x))

ids <- data.frame("FID" = df$ID, "IID" = df$ID)
df <- cbind(ids, df)

#Write tables of samples to keep: 
subset_df <- df[df$check_se %in% c("brit", "white", "recode.white"), c("FID", "IID")]
write.table(subset_df, "../sex.het.enrichment/all.white.comb.samples.to.keep", quote = F, row.names = F, col.names = F, sep = " ")

subset_df <- df[df$check_se %in% c("brit", "white", "recode.white") & df$Submitted_Gender == "M", c("FID", "IID")]
write.table(subset_df, "../sex.het.enrichment/all.white.men.samples.to.keep", quote = F, row.names = F, col.names = F, sep = " ")

subset_df <- df[df$check_se %in% c("brit", "white", "recode.white") & df$Submitted_Gender == "F", c("FID", "IID")]
write.table(subset_df, "../sex.het.enrichment/all.white.women.samples.to.keep", quote = F, row.names = F, col.names = F, sep = " ")

#Write a phenotype file in PLINK format
write.table(df, "../sex.het.enrichment/plink.phenotype.file.from.ukbb.cases.for.logistic.180921.txt", 
	quote = F, row.names = F, sep = " ")

