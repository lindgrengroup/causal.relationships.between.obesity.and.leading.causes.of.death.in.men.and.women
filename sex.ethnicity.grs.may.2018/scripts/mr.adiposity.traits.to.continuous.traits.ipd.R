#!/bin/env Rscript
#$-cwd

library(mada)

qc.pregnant <- function(data) {
	cleaned <- subset(data, is.na(data$pregnant) | data$pregnant == 0)
	return(cleaned)
}

qc.below.15.bmi <- function(data) {
	cleaned <- subset(data, is.na(data$bmi) | data$bmi > 15)
	return(cleaned)
}

run.regression.wald <- function(trait_df, trait, exposure_name, outcome_name, extra_adjustment) {
	pcs <- paste("PC", 1:10, sep = "", collapse = " + ")
        if (exposure_name %in% c("res_whr_inv", "exposure_sd") & extra_adjustment == "-") {
                linreg.form <- as.formula(paste(exposure_name, "~",
                                "SCORESUM+dummy_array+", pcs, sep = " "))
        } else if (extra_adjustment == "-") {
                linreg.form <- as.formula(paste(exposure_name, "~",
                                "SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                pcs, sep = " "))
        } else if (exposure_name %in% c("res_whr_inv", "exposure_sd") & extra_adjustment == "smoking") {
                linreg.form <- as.formula(paste(exposure_name, "~",
                                "SCORESUM+smoker_cases+dummy_array+", pcs, sep = " "))
	} else if (extra_adjustment == "smoking") {
		linreg.form <- as.formula(paste(exposure_name, "~",
                                "SCORESUM+smoker_cases+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                pcs, sep = " "))
	}
	print(linreg.form)
	raw_fit1 <- lm(linreg.form, data = trait_df, na.action = na.exclude)
        if (outcome_name %in% c("res_whr_inv", "outcome_sd") & extra_adjustment == "-") {
                linreg.form <- as.formula(paste(outcome_name, "~",
                                "SCORESUM+dummy_array+", pcs, sep = " "))
        } else if (extra_adjustment == "-") {
                linreg.form <- as.formula(paste(outcome_name, "~",
                                "SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                pcs, sep = " "))
        } else if (outcome_name %in% c("res_whr_inv", "outcome_sd") & extra_adjustment == "smoking") {
		linreg.form <- as.formula(paste(outcome_name, "~",
                                "SCORESUM+smoker_cases+dummy_array+", pcs, sep = " "))
	} else if (extra_adjustment == "smoking") {
		linreg.form <- as.formula(paste(outcome_name, "~",
                                "SCORESUM+smoker_cases+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                pcs, sep = " "))
	}
	print(linreg.form)
	raw_fit2 <- lm(linreg.form, data = trait_df, na.action = na.exclude)

	#Compute Beta and SEs
	ratio_beta <- summary(raw_fit2)$coef["SCORESUM", "Estimate"]/summary(raw_fit1)$coef["SCORESUM", "Estimate"]
	ratio_se <- sqrt( (summary(raw_fit2)$coef["SCORESUM", "Std. Error"]^2/summary(raw_fit1)$coef["SCORESUM", "Estimate"]^2) +
                        ((summary(raw_fit2)$coef["SCORESUM", "Estimate"]^2 * summary(raw_fit1)$coef["SCORESUM", "Std. Error"]^2) /
                                summary(raw_fit1)$coef["SCORESUM", "Estimate"]^4))
			
        ratio_z <- ratio_beta/ratio_se
        ratio_p <- 2*pnorm(abs(ratio_z), lower.tail = F)

        fit_raw <- data.frame("beta" = ratio_beta, "se" = ratio_se, "p" = ratio_p)
        return(fit_raw)
}

current_row <- function(results_df, snp_group, grs_unit, eth_group, sex_group, grs_trait, trait, file,
		fit_raw, exposure_unit, outcome_unit, trait_df, function_name, extra_adjustment) {

	current_row <- data.frame("snp_group" = snp_group, "grs_unit" = grs_unit, 
        		"eth_group" = eth_group, "sex_group" = sex_group,
			"exposure_unit" = exposure_unit, "outcome_unit" = outcome_unit,
                        "trait" = ifelse(extra_adjustment == "smoking", paste(trait, "_smoking", sep = ""), trait), 
			"grs_trait" = grs_trait,
                        "grs_p" = fit_raw[, "p"],
                        "grs_beta" = fit_raw[, "beta"],
                        "grs_se" = round(fit_raw[, "se"], 4),
                        "grs_lci_beta" = round(fit_raw[, "beta"] - (1.96 * fit_raw[ , "se"]), 4),
                        "grs_uci_beta" = round(fit_raw[ , "beta"] + (1.96 * fit_raw[ , "se"]), 4),
                        "sd_exposure_corresponds_to" = sd(trait_df[[grs_trait]], na.rm = T),
			"sd_outcome_corresponds_to" = sd(trait_df[[trait]], na.rm = T),
                        "grs_file" = file,
                        "n_complete_cases" = NA,
			"function_name" = function_name,
			"extra_adjustment" = extra_adjustment,
                        stringsAsFactors = F)

	results_df <- merge(results_df, current_row, by.x = colnames(results_df),
                            by.y = colnames(current_row), all = T)
	return(results_df)	
}

tsls.cont.mr.grs <- function(pheno, grs_results, snp_groups, eth_groups, grs_unit_groups, sex_groups, traits) {

	results_df <-  data.frame("snp_group" = character(), "grs_unit" = character(), 
                        "eth_group" = character(), "sex_group" = character(),
                        "exposure_unit" = character(), "outcome_unit" = character(),
                        "trait" = character(), "grs_trait" = character(),
                        "grs_p" = numeric(), "grs_beta" = numeric(),
                        "grs_se" = numeric(), "grs_lci_beta" = numeric(),
                        "grs_uci_beta" = numeric(), 
                        "sd_exposure_corresponds_to" = numeric(),
                        "sd_outcome_corresponds_to" = numeric(),
                        "grs_file" = character(), 
                        "n_complete_cases" = integer(),
			"function_name" = character(),
			"extra_adjustment" = character(), 
                        stringsAsFactors = F)

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

                                #Make columns with WHR adjusted for BMI
                                whrfit <- lm(whr ~ age_assessment + age_squared + dummy_sex + bmi + assessment_centre,
                                        data = sex_df, na.action = na.exclude)
                                sex_df$res_whr <- residuals(whrfit)
                                sex_df$res_whr_inv <- qnorm((rank(sex_df$res_whr,
                                        na.last = "keep") - 0.5) / sum(!is.na(sex_df$res_whr)))

                                for (grs_unit in grs_unit_groups) {
                                        grs_unit_df <- sex_df
                                        if (grs_unit == "sd_scoresum") {
                                                grs_unit_df$SCORESUM <- (grs_unit_df$SCORESUM - mean(grs_unit_df$SCORESUM,
                                                                        na.rm = T))/sd(grs_unit_df$SCORESUM, na.rm = T)
                                        }

                                        #Make the formula:
                                        for (trait in traits) {
						for (extra_adjustment in c("-", "smoking")) {
                                                        grs_trait <- gsub("\\.(.)+", "", snp_group)

                                                        #Decide whether to run based on the sig analysis
							decision_sex_group <- gsub("(.)+.eur\\.", "", snp_group)
                                                        decision_sex_group <- gsub("\\.(.)+", "", decision_sex_group)
                                                        decision_snp_group <- paste(grs_trait, ".eur.", decision_sex_group, ".pulit.sig", sep = "")

							if (grs_trait == "whradjbmi") {
								grs_trait <- "res_whr_inv"
							}

							#Decide which traits to run based on if any sex significant in all.white or brit.irish in the GRS analysis
							if (grs_trait == trait | nrow(grs_results[grs_results$snp_group == decision_snp_group & 
								grs_results$eth_group %in% c("all.white", "brit.irish") &
								grs_results$sex_group == sex_group &
								grs_results$grs_unit == grs_unit &
								grs_results$trait == trait, ]) < 1 ) {
								next()
							}

                                                        #Make RINTed traits - if WHRadjBMI, use that one directly
							trait_df <- grs_unit_df
                                                        if (grs_trait != "res_whr_inv") {
                                                                rint.formula.exp <- as.formula(paste(grs_trait,
                                                                "~age_assessment+age_squared+dummy_sex+assessment_centre", sep = " "))
                                                                rint_fit_exp <- lm(rint.formula.exp, data = trait_df, na.action = na.exclude)
                                                                trait_df$res_rint_exp <- residuals(rint_fit_exp)
                                                                trait_df$exposure_sd <- qnorm((rank(trait_df$res_rint_exp,
                                                                                na.last = "keep") - 0.5) / sum(!is.na(trait_df$res_rint_exp)))
                                                        } else if (grs_trait == "res_whr_inv") {
                                                                trait_df$exposure_sd <- trait_df$res_whr_inv
                                                        }
							print(grs_trait)
							print(rint_fit_exp)
							if (trait != "res_whr_inv") {
                                                                rint.formula.out <- as.formula(paste(trait,
                                                                        "~age_assessment+age_squared+dummy_sex+assessment_centre", sep = " "))
                                                                        rint_fit_out <- lm(rint.formula.out, data = trait_df, na.action = na.exclude)
                                                                        trait_df$res_rint_out <- residuals(rint_fit_out)
                                                                        trait_df$outcome_sd <- qnorm((rank(trait_df$res_rint_out,
                                                                                na.last = "keep") - 0.5) / sum(!is.na(trait_df$res_rint_out)))
                                                        } else if (trait == "res_whr_inv") {
                                                                trait_df$outcome_sd <- trait_df$res_whr_inv
                                                        }
							print(trait)
							print(rint_fit_out)
								
							#For GRS-trait in clinical units, outcome in clinical units
							fit_raw <- run.regression.wald(trait_df, trait, exposure_name = grs_trait, outcome_name = trait, 
									extra_adjustment = extra_adjustment)
							results_df <- current_row(results_df, snp_group, grs_unit, eth_group, sex_group, 
									grs_trait, trait, file,
									fit_raw, exposure_unit = "clin", outcome_unit = "clin", trait_df, "wald", extra_adjustment = extra_adjustment)

							#For GRS-trait in SD-units, outcome in clinical units
							fit_raw <- run.regression.wald(trait_df, trait, exposure_name = "exposure_sd", outcome_name = trait, 
								extra_adjustment = extra_adjustment)
							results_df <- current_row(results_df, snp_group, grs_unit, eth_group, sex_group,
                                                       	        grs_trait, trait, file,
                                                               	fit_raw, exposure_unit = "sd", outcome_unit = "clin", trait_df, "wald", extra_adjustment = extra_adjustment)
						
							#For GRS-trait in clinical units, outcome in SD-units
							fit_raw <- run.regression.wald(trait_df, trait, exposure_name = grs_trait, outcome_name = "outcome_sd", 
								extra_adjustment = extra_adjustment)
                	                       	        results_df <- current_row(results_df, snp_group, grs_unit, eth_group, sex_group, 
								grs_trait, trait, file,
                                	                               fit_raw, exposure_unit = "clin", outcome_unit = "sd", trait_df, "wald", extra_adjustment = extra_adjustment)

							#For GRS-trait in SD-units, outcome in SD-units
							fit_raw <- run.regression.wald(trait_df, trait, exposure_name = "exposure_sd", outcome_name = "outcome_sd", 
									extra_adjustment = extra_adjustment)
                                                	results_df <- current_row(results_df, snp_group, grs_unit, eth_group, sex_group,
       	                                                        grs_trait, trait, file,
               	                                                fit_raw, exposure_unit = "sd", outcome_unit = "sd", trait_df, "wald", extra_adjustment = extra_adjustment)

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
pheno <- read.table("/well/lindgren/jc/ukbb/ukbb.samples.passing.qc.relevant.pheno.plus.diagnoses.180509.txt",
                 stringsAsFactors = F, header = T)
pheno$age_squared <- pheno$age_assessment * pheno$age_assessment
pheno$assessment_centre <- factor(pheno$assessment_centre)

smoking <- read.table("../ukbb_cases_for_logistic_regression_180620.txt", stringsAsFactors = F, header = T, sep = " ")
smoking <- smoking[, c("ID", "smoker_cases")]
smoking$smoker_cases <- factor(smoking$smoker_cases)

pheno <- merge(pheno, smoking, by.x = "ID", by.y = "ID")

#Make blood pressure columns - mean of baseline
pheno$sbp <- apply(pheno[, c("X4080.0_0", "X4080.0_1")], 1, function(x) {mean(x, na.rm = T)})
pheno$dbp <- apply(pheno[, c("X4079.0_0", "X4079.0_1")], 1, function(x) {mean(x, na.rm = T)})

#Increase the blood pressure if on blood-pressure lowering medications at baseline:
pheno[apply(pheno[, grep("X6177.0_|X6153.0_", names(pheno))], 1,
        function(x) {any(x == 2, na.rm = T)}), "sbp"] <- pheno[apply(pheno[, grep("X6177.0_|X6153.0_",
        names(pheno))], 1, function(x) {any(x == 2, na.rm = T)}), "sbp"] + 15

pheno[apply(pheno[, grep("X6177.0_|X6153.0_", names(pheno))], 1,
        function(x) {any(x == 2, na.rm = T)}), "dbp"] <- pheno[apply(pheno[, grep("X6177.0_|X6153.0_",
        names(pheno))], 1, function(x) {any(x == 2, na.rm = T)}), "dbp"] + 10

pheno$dummy_array <- ifelse(pheno$genotyping_array == "UKBB", 0, 1)
pheno$dummy_sex <- ifelse(pheno$Submitted_Gender == "F", 0, 1)

#Clean the phenotype file:
pheno <- qc.pregnant(pheno)
pheno <- qc.below.15.bmi(pheno)

#Get the results from the GRSs that are significicant:
grs_results <- read.table("../results.linear.regressions.180514/anthro.results.table.180521.txt", 
		stringsAsFactors = F, header = T, sep = "\t")

#Remove the adjusted for smoking analyses
adjusted_smoking <- unique(grs_results$trait[grep("_smoking", grs_results$trait)])
grs_results <- grs_results[!(grs_results$trait %in% adjusted_smoking), ]

comb_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", "whradjbmi.eur.comb.pulit.sig")
male_groups <- c("bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig")
female_groups <- c("bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.sig")

#Subset to relevant analyses
grs_results <- grs_results[((grs_results$snp_group %in% c(comb_groups) &
		grs_results$sex_group == "comb") |
		(grs_results$snp_group %in% c(male_groups) &
		grs_results$sex_group == "men") |
		(grs_results$snp_group %in% c(female_groups) &
                grs_results$sex_group == "women")) &
		grs_results$grs_unit == "raw_scoresum" &
		grs_results$eth_group %in% c("all.white", "brit.irish"), ]

#Get the significant P-value, Bonferroni-adjustment for GRS traits and outcome traits
sig_p <- 0.05/(length(unique(grs_results$trait)) * length(unique(gsub("\\..*", "", grs_results$snp_group))))

#Remove the comb regressions that aren't significant
for (i in 1:length(comb_groups)) {
        for (grs_unit in unique(grs_results$grs_unit)) {
                for (trait in unique(grs_results$trait)) {
                         if (nrow(grs_results[(grs_results$snp_group %in% comb_groups[i]) &
                                  grs_results$grs_unit == grs_unit & grs_results$trait == trait & grs_results$grs_p < sig_p, ]) == 0) {

                                grs_results <- grs_results[!(grs_results$snp_group %in% c(comb_groups[i]) &
                                                grs_results$grs_unit == grs_unit & grs_results$trait == trait), ]
                        }
                }
        }
}

#Remove the sex-specific regression where neither sex is significant
for (i in 1:length(female_groups)) {
	for (grs_unit in unique(grs_results$grs_unit)) {
		for (trait in unique(grs_results$trait)) {
			if (nrow(grs_results[grs_results$snp_group %in% c(female_groups[i], male_groups[i]) &
				grs_results$grs_unit == grs_unit & grs_results$trait == trait & grs_results$grs_p < sig_p, ]) == 0) {
					
						grs_results <- grs_results[!(grs_results$snp_group %in% c(female_groups[i], male_groups[i]) &
                               		       grs_results$grs_unit == grs_unit & grs_results$trait == trait), ]
			}
		}
	}
}

#Create the vectors to loop over
snp_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig",
                "whradjbmi.eur.comb.pulit.sig",
                "bmi.eur.men.pulit", "bmi.eur.men.pulit.sig", "bmi.eur.men.pulit.phet",
                "whr.eur.men.pulit", "whr.eur.men.pulit.sig", "whr.eur.men.pulit.phet",
                "whradjbmi.eur.men.pulit", "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.phet",
                "bmi.eur.women.pulit", "bmi.eur.women.pulit.sig", "bmi.eur.women.pulit.phet",
                "whr.eur.women.pulit", "whr.eur.women.pulit.sig", "whr.eur.women.pulit.phet",
                "whradjbmi.eur.women.pulit", "whradjbmi.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.phet",
                "bmi.eur.men.0.01.fdr", "bmi.eur.men.0.05.fdr",
                "bmi.eur.men.0.1.fdr", "bmi.eur.women.0.01.fdr",
                "bmi.eur.women.0.05.fdr", "bmi.eur.women.0.1.fdr",
                "whr.eur.men.0.01.fdr", "whr.eur.men.0.05.fdr",
                "whr.eur.men.0.1.fdr", "whr.eur.women.0.01.fdr",
                "whr.eur.women.0.05.fdr", "whr.eur.women.0.1.fdr",
                "whradjbmi.eur.men.0.01.fdr", "whradjbmi.eur.men.0.05.fdr",
                "whradjbmi.eur.men.0.1.fdr", "whradjbmi.eur.women.0.01.fdr",
                "whradjbmi.eur.women.0.05.fdr", "whradjbmi.eur.women.0.1.fdr",
                "bmi.eur.comb.giukbb.sig", "whr.eur.comb.giukbb.sig", "whradjbmi.eur.comb.giukbb.sig",
                "bmi.eur.men.giukbb.sig", "whr.eur.men.giukbb.sig", "whradjbmi.eur.men.giukbb.sig",
                "bmi.eur.women.giukbb.sig", "whr.eur.women.giukbb.sig", "whradjbmi.eur.women.giukbb.sig", 
                "bmi.eur.comb.internal.unweighted", "bmi.eur.women.internal.unweighted", "bmi.eur.men.internal.unweighted",
                "whr.eur.comb.internal.unweighted", "whr.eur.women.internal.unweighted", "whr.eur.men.internal.unweighted",
                "whradjbmi.eur.comb.internal.unweighted", "whradjbmi.eur.women.internal.unweighted", "whradjbmi.eur.men.internal.unweighted")
eth_groups <- c("all.white", "brit.irish")

#Get same results for raw scoresum and sd_scoresum
grs_unit_groups <- c("raw_scoresum") 
sex_groups <- c("comb", "women", "men")
traits <- c("sbp", "dbp", "whr", "wc", "hc", "bmi", "res_whr_inv")

results_df <- tsls.cont.mr.grs(pheno, grs_results, snp_groups, eth_groups, grs_unit_groups, sex_groups, traits)

write.table(results_df, "../results.mr.180730/ipd.mr.continuous.results.180815.EXTRA.txt", sep = "\t", row.names = F, quote = F)

sink("../results.mr.180730/ipd.mr.continuous.results.180815.WARNINGS.txt", append = F)
print(warnings())
sink()

#Check for heterogeneity
results_df$cochran_weights <- 1 / (results_df$grs_se^2)

women_snp_groups <- snp_groups[grep("\\.women", snp_groups)]
men_snp_groups <- snp_groups[grep("\\.men", snp_groups)]

women_snp_groups <- women_snp_groups[order(women_snp_groups)]
men_snp_groups <- men_snp_groups[order(men_snp_groups)]

for (i in 1:length(women_snp_groups)) {
        for (eth_group in unique(results_df$eth_group)) {
                for (trait in unique(results_df$trait)) {
                        for (exposure_unit in unique(results_df$exposure_unit)) {
				for (outcome_unit in unique(results_df$outcome_unit)) {

					general <- results_df$eth_group == eth_group &
                                                results_df$trait == trait &
                                                results_df$exposure_unit == exposure_unit &
						results_df$outcome_unit == outcome_unit

					if (nrow(results_df[general & results_df$snp_group %in% c(women_snp_groups[i], men_snp_groups[i]) &
                                                results_df$sex_group %in% c("women", "men"), ]) == 2) {
		
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
}

results_df$cochrans_i2 <- (results_df$cochrans_q - 1)/results_df$cochrans_q
results_df$cochrans_i2[results_df$cochrans_i2 < 0] <- 0

write.table(results_df, "../results.mr.180730/ipd.mr.continuous.results.180815.txt", sep = "\t", row.names = F, quote = F)




