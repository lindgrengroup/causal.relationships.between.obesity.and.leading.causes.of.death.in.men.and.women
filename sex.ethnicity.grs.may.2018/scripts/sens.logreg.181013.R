#!/bin/env Rscript
#$ -cwd

##########################################################
###Sex and ethnicity WHR GRS on anthropometric traits ####
############ Jenny Censin, 2018-05-20 ###################

library(mada)

logreg <- function(pheno, snp_groups, eth_groups, sex_groups, case_columns) {

	results_df <- data.frame("snp_group" = character(), "grs_unit" = character(), "eth_group" = character(), 
				"sex_group" = character(), "case_column" = character(), 
                                "grs_p" = numeric(), "grs_beta" = numeric(),
				"grs_se" = numeric(), "grs_lci_beta" = numeric(), "grs_uci_beta" = numeric(),
				"grs_or" = numeric(), "grs_lci_or" = numeric(), "grs_uci_or" = numeric(),
				"grs_file" = character(), "n_cases" = integer(), "n_controls" = integer(),
				"n_excluded" = integer(), 
				"n_complete_obs" = integer(), 
				"extra_adjustment" = character(), 
				stringsAsFactors = F)

	#Make the covariates
	pcs <- paste("PC", 1:10, sep = "", collapse = " + ")

	for (snp_group in snp_groups) {
	        file <- paste("../bmi.whr.profiles/grs.", snp_group, ".180513.profile", sep = "")
        	scores <- read.table(file, stringsAsFactors = F, header = T)
		df <- merge(x = scores, y = pheno, by.x = "IID", by.y = "ID")
		
        	for (eth_group in eth_groups) {
                	if (eth_group == "brit.irish") {
                        	eth_df <- subset(df, df$check_se %in% c("brit"))
                	} else if (eth_group == "all.white") {
                        	eth_df <- subset(df, df$check_se %in% c("brit", "white", "recode.white"))
                	}
              		for (sex_group in sex_groups) {
                        	if (sex_group == "comb") {
                                	sex_df <- eth_df
                        	} else if (sex_group == "women") {
                                	sex_df <- subset(eth_df, eth_df$Submitted_Gender == "F")
                        	} else if (sex_group == "men") {
                                	sex_df <- subset(eth_df, eth_df$Submitted_Gender == "M")
                        	}

				for (grs_unit in c("raw_scoresum")) {
                                        grs_unit_df <- sex_df
                                        if (grs_unit == "sd_scoresum") {
                                                grs_unit_df$SCORESUM <- (grs_unit_df$SCORESUM - mean(grs_unit_df$SCORESUM,
                                                                        na.rm = T))/sd(grs_unit_df$SCORESUM, na.rm = T)
                                        }

					for (case_column in case_columns) {
						for (extra_adjustment in c("-")) {
							if (extra_adjustment == "smoking" & case_column == "smoker_cases") {
                	                                                next()
							}
	
							#Exclude the cases that are to be excluded
							n_excluded <- nrow(grs_unit_df[grs_unit_df[, case_column] == 2, ])
							case_df <- grs_unit_df[grs_unit_df[, case_column] != 2, ]
							case_df[, case_column] <- factor(case_df[, case_column])
						
							if (nrow(case_df) < 1) {							
								next
							}
	
                                                        grs_sex <- gsub("(.)+\\.eur\\.", "", snp_group)
                                                        grs_sex <- gsub("\\.(.)*", "", grs_sex)

                                                        if (grs_sex != "comb" & grs_sex != sex_group) {
                                                                next()
                                                        }
								
							#Make the formula:
							if (extra_adjustment == "-") {
	        		                                logreg.form <- as.formula(paste(case_column, 
									"~SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                		                                        pcs, sep = " "))
							} else if (extra_adjustment == "smoking") {
								case_df$smoker_cases <- as.factor(case_df$smoker_cases)
								logreg.form <- as.formula(paste(case_column, 
									"~SCORESUM+smoker_cases+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                                        	        pcs, sep = " "))
							}

        	                	                fit_raw <- glm(formula = logreg.form, data = case_df, na.action = na.exclude, family = binomial(link="logit"))

							print(summary(fit_raw))
                        	                	fit <- as.data.frame(summary(fit_raw)$coefficients)	
			
							current_row <- data.frame("snp_group" = snp_group, "grs_unit" = grs_unit, 
								"eth_group" = eth_group, "sex_group" = sex_group, 
								"case_column" = ifelse(extra_adjustment == "smoking", paste(case_column, "_smoking", sep = ""), case_column),
								"grs_p" = fit["SCORESUM", "Pr(>|z|)"], 
								"grs_beta" = fit["SCORESUM", "Estimate"],
								"grs_se" = round(fit["SCORESUM", "Std. Error"], 4),
								"grs_lci_beta" = round(fit["SCORESUM", "Estimate"] - (1.96 * fit["SCORESUM", "Std. Error"]), 4),
								"grs_uci_beta" = round(fit["SCORESUM", "Estimate"] + (1.96 * fit["SCORESUM", "Std. Error"]), 4), 
								"grs_or" = round(exp(fit["SCORESUM", "Estimate"]), 4),
								"grs_lci_or" = round(exp(fit["SCORESUM", "Estimate"] - (1.96 * fit["SCORESUM", "Std. Error"])), 4), 
								"grs_uci_or" = round(exp(fit["SCORESUM", "Estimate"] + (1.96 * fit["SCORESUM", "Std. Error"])), 4),
								"grs_file" = file, 
								"n_cases" = nrow(case_df[case_df[case_column] == 1, ]),
								"n_controls" = nrow(case_df[case_df[case_column] == 0, ]), 
								"n_excluded" = n_excluded,
								"n_complete_obs" = nobs(fit_raw),
								"extra_adjustment" = extra_adjustment,
								stringsAsFactors = F)
							print(current_row)
	
							results_df <- merge(results_df, current_row, by.x = colnames(results_df), 
								by.y = colnames(current_row), all = T)
							print(results_df)
						}	
					}
				}
			}
		}
	}
	return(results_df)	
}

#Read in phenotype file:
pheno <- read.table("../ukbb_cases_for_logistic_regression_180620.txt",
                 stringsAsFactors = F, header = T, sep = " ")
pheno$age_squared <- pheno$age_assessment * pheno$age_assessment
pheno$dummy_array <- ifelse(pheno$genotyping_array == "UKBB", 0, 1)
pheno$dummy_sex <- ifelse(pheno$Submitted_Gender == "F", 0, 1)
pheno$assessment_centre <- factor(pheno$assessment_centre)

#Create the vectors to loop over
eth_groups <- c("all.white")
sex_groups <- c("comb", "women", "men")
case_columns <- c("t2d_cases_probposs")
types_internal <- c("normal", "half", "randomnoise", "systematicnoise", "unweighted")
types_fdr <- c("0.01.fdr", "0.05.fdr", "0.1.fdr")
x <- paste(rep(c("bmi", "whr", "whradjbmi"), each = length(sex_groups)), "eur", sex_groups, sep = ".")
internal <- paste(rep(x, each = length(types_internal)), "internal", types_internal, sep = ".")
y <- paste(rep(c("bmi", "whr", "whradjbmi"), each = length(c("men", "women"))), "eur", c("men", "women"), sep = ".")
fdr <- paste(rep(y, each = length(types_fdr)), types_fdr, sep = ".")
snp_groups <- c(internal, fdr)

results_df <- logreg(pheno, snp_groups, eth_groups, sex_groups, case_columns)

write.table(results_df, "../results.logistic.regressions.180514/sens.log.results.table.180627.EXTRA.txt",
                row.names = F, sep = "\t", quote = F)

sink("../results.logistic.regressions.180514/sens.log.results.table.180627.WARNINGS.txt", append = F)
print(warnings())
sink()

comb_snp_groups <- snp_groups[grep("\\.comb\\.", snp_groups)]
women_snp_groups <- snp_groups[grep("\\.women\\.", snp_groups)]
men_snp_groups <- snp_groups[grep("\\.men\\.", snp_groups)]

women_snp_groups <- women_snp_groups[order(women_snp_groups)]
men_snp_groups <- men_snp_groups[order(men_snp_groups)]

#Prepare some columns required for Heterogeneity calculations
results_df$cochran_weights <- 1 / (results_df$grs_se^2)

#Adjustment for smoking is taken cared of by the fact that I've amended with "_smoking" in those cases
for (i in 1:length(women_snp_groups)) {
        for (eth_group in unique(results_df$eth_group)) {
                for (case_column in unique(results_df$case_column)) {
                        for (grs_unit in unique(results_df$grs_unit)) {
                                        general <- results_df$eth_group == eth_group &
                                                results_df$case_column == case_column &
                                                results_df$grs_unit == grs_unit

					  if (i <= length(comb_snp_groups) &
                                                nrow(results_df[general & results_df$snp_group == comb_snp_groups[i] &
                                                results_df$sex_group %in% c("women", "men"), ]) == 2) {

                                                results_df$cochrans_q[(general & results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group %in% c("women", "men"))] <- as.numeric(cochran.Q(c(results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "women",
                                                        "grs_beta"],
                                                        results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "men",
                                                        "grs_beta"]),
                                                        c(results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "women",
                                                        "cochran_weights"],
                                                        results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "men",
                                                        "cochran_weights"]))[1])

                                                results_df$cochrans_p[(general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group %in% c("women", "men"))] <- as.numeric(cochran.Q(c(results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "women",
                                                        "grs_beta"],
                                                        results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "men",
                                                        "grs_beta"]),
                                                        c(results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "women",
                                                        "cochran_weights"],
                                                        results_df[general &
                                                        results_df$snp_group == comb_snp_groups[i] &
                                                        results_df$sex_group == "men",
                                                        "cochran_weights"]))[2])

                                        }

                                if (nrow(results_df[(general &
                                        results_df$snp_group == women_snp_groups[i] &
                                        results_df$sex_group == "women") |
                                        (general &
                                        results_df$snp_group == men_snp_groups[i] &
                                        results_df$sex_group == "men"), ]) == 2) {


                                        results_df$cochrans_q[(general &
                                                results_df$snp_group == women_snp_groups[i] &
                                                results_df$sex_group == "women") |
                                                (general &
                                                results_df$snp_group == men_snp_groups[i] &
                                                results_df$sex_group == "men")] <- as.numeric(cochran.Q(c(results_df[general &
                                                results_df$snp_group == women_snp_groups[i] &
                                                results_df$sex_group == "women",
                                                "grs_beta"],
                                                results_df[general &
                                                results_df$snp_group == men_snp_groups[i] &
                                                results_df$sex_group == "men",
                                                "grs_beta"]),
                                                c(results_df[general &
                                                results_df$snp_group == women_snp_groups[i] &
                                                results_df$sex_group == "women",
                                                "cochran_weights"],
                                                results_df[general &
                                                results_df$snp_group == men_snp_groups[i] &
                                                results_df$sex_group == "men",
                                                "cochran_weights"]))[1])

                                        results_df$cochrans_p[(general &
                                                results_df$snp_group == women_snp_groups[i] &
                                                results_df$sex_group == "women") |
                                                (general &
                                                results_df$snp_group == men_snp_groups[i] &
                                                results_df$sex_group == "men")] <- as.numeric(cochran.Q(c(results_df[general &
                                                results_df$snp_group == women_snp_groups[i] &
                                                results_df$sex_group == "women",
                                                "grs_beta"],
                                                results_df[general &
                                                results_df$snp_group == men_snp_groups[i] &
                                                results_df$sex_group == "men",
                                                "grs_beta"]),
                                                c(results_df[general &
                                                results_df$snp_group == women_snp_groups[i] &
                                                results_df$sex_group == "women",
                                                "cochran_weights"],
                                                results_df[general &
                                                results_df$snp_group == men_snp_groups[i] &
                                                results_df$sex_group == "men",
                                                "cochran_weights"]))[2])

				}                                
                        }
                }
        }
}

results_df$cochrans_i2 <- (results_df$cochrans_q - 1)/results_df$cochrans_q

results_df$cochrans_i2[results_df$cochrans_i2 < 0] <- 0

write.table(results_df, "../results.logistic.regressions.180514/sens.log.results.table.180627.txt", 
		row.names = F, sep = "\t", quote = F)




