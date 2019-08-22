#!/bin/bash
#$-cwd

#Quality control functions
qc.pregnant <- function(data) {
cleaned <- subset(data, is.na(data$pregnant) | data$pregnant == 0)
return(cleaned)
}

qc.below.15.bmi <- function(data) {
cleaned <- subset(data, is.na(data$bmi) | data$bmi > 15)
return(cleaned)
}

#Read in phenotype file
df <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
                 stringsAsFactors = F, header = T)
df <- qc.pregnant(df)
df <- qc.below.15.bmi(df)

#Make PCs
pcs <- paste("PC", 1:10, sep = "", collapse = " + ")

#Subset to Europeans only and create dummy variables
df <- subset(df, df$check_se %in% c("brit", "white", "recode.white"))
df$dummy_array <- ifelse(df$genotyping_array == "UKBB", 0, 1)
df$dummy_sex <- ifelse(df$Submitted_Gender == "F", 0, 1)
colnames(df)[colnames(df) == "ID"] <- "IID" 
df$FID <- df$IID
df$assessment_centre <- factor(df$assessment_centre)

#Create vectors to loop over
sex_groups <- c("comb", "women", "men")
traits <- c("bmi", "whr", "res_whr_inv")

#Get residuals
for (sex_group in sex_groups) {
	if (sex_group == "comb") {
        	sex_df <- df
        } else if (sex_group == "women") {
                sex_df <- df[df$Submitted_Gender == "F", ]
        } else if (sex_group == "men") {
                sex_df <- df[df$Submitted_Gender == "M", ]
        }
	for (trait in traits) {
		if (trait == "res_whr_inv") {
			linreg.form <- as.formula(paste("whr ~ age_assessment + age_squared + dummy_sex + dummy_array + bmi + assessment_centre +", pcs, sep = ""))
		} else {
			linreg.form <- as.formula(paste(trait, "~ age_assessment + age_squared + dummy_sex + dummy_array + assessment_centre +", pcs, sep = ""))
		}
			
		fit <- lm(formula = linreg.form, data = sex_df, na.action = na.exclude)
 		print(fit)
		column_name <- paste(trait, "_", sex_group, sep = "")
       		print(column_name)
		sex_df[, column_name] <- residuals(fit)
		sex_df[, column_name] <- qnorm((rank(sex_df[, column_name],
                                        na.last = "keep") - 0.5) / sum(!is.na(sex_df[, column_name])))

	}
	sex_df <- sex_df[, grep(paste("^IID$|_", sex_group, "$", sep = ""), colnames(sex_df))]
	df <- merge(df, sex_df, by.x = "IID", by.y = "IID", all.x = T)
}

df <- df[, grep(paste("^IID$|^FID$|^Submitted_Gender$|_comb$|_women$|_men$", sep = ""), colnames(df))]
df <- df[, c("FID", "IID", "Submitted_Gender", colnames(df)[colnames(df) %in% colnames(df)[grep("_comb$|_women$|_men$", colnames(df))]])]

write.table(df, "../ukbb.residuals.for.adiposity.traits.180906.txt", quote = F, sep = "\t", row.names = F)

