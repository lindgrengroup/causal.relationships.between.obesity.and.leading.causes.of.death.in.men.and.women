#!/bin/env Rscript
#$-cwd

library(mada)
library(ivmodel)

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
        if (exposure_name %in% c("res_whr_inv", "exposure_sd") & extra_adjustment == "smoking") {
                linreg.form <- as.formula(paste(exposure_name, "~",
                                "SCORESUM+smoker_cases+dummy_array+", pcs, sep = " "))
        } else if (exposure_name %in% c("res_whr_inv", "exposure_sd") & extra_adjustment == "-") {
		linreg.form <- as.formula(paste(exposure_name, "~",
                                "SCORESUM+dummy_array+", pcs, sep = " "))
	} 
        raw_fit1 <- lm(linreg.form, data = trait_df[trait_df[[outcome_name]] != 1, ], na.action = na.exclude)

	if (extra_adjustment == "smoking") {
	        linreg.form <- as.formula(paste(outcome_name, "~",
                                "SCORESUM+smoker_cases+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                pcs, sep = " "))
	} else if (extra_adjustment == "-") {
		linreg.form <- as.formula(paste(outcome_name, "~",
                                "SCORESUM+dummy_array+dummy_sex+age_assessment+age_squared+assessment_centre+",
                                pcs, sep = " "))
	}

        raw_fit2 <- glm(linreg.form, data = trait_df, na.action = na.exclude, family = binomial(link="logit"))

	#Two first terms of delta method
	ratio_beta <- summary(raw_fit2)$coef["SCORESUM", "Estimate"]/summary(raw_fit1)$coef["SCORESUM", "Estimate"]
	ratio_se <- sqrt( (summary(raw_fit2)$coef["SCORESUM", "Std. Error"]^2/summary(raw_fit1)$coef["SCORESUM", "Estimate"]^2) +
			((summary(raw_fit2)$coef["SCORESUM", "Estimate"]^2 * summary(raw_fit1)$coef["SCORESUM", "Std. Error"]^2) /
				summary(raw_fit1)$coef["SCORESUM", "Estimate"]^4))
	
	ratio_z <- ratio_beta/ratio_se
	ratio_p <- 2*pnorm(abs(ratio_z), lower.tail = F)

	fit_raw <- data.frame("beta" = ratio_beta, "se" = ratio_se, "p" = ratio_p)
        return(fit_raw)
}

current_row_wald <- function(results_df, snp_group, grs_unit, eth_group, sex_group, grs_trait, trait, file,
                fit_raw, exposure_unit, outcome_unit, trait_df, function_name, n_excluded, extra_adjustment) {

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
                        "grs_or" = round(exp(fit_raw[, "beta"]), 4),
                        "grs_lci_or" = round(exp(fit_raw[, "beta"] - (1.96 * fit_raw[, "se"])), 4),
                        "grs_uci_or" = round(exp(fit_raw[, "beta"] + (1.96 * fit_raw[, "se"])), 4),
                        "n_cases" = nrow(trait_df[trait_df[trait] == 1, ]),
                        "n_controls" = nrow(trait_df[trait_df[trait] == 0, ]),
                        "n_excluded" = n_excluded,
                        "grs_file" = file,
                        "n_complete_cases" = nrow(trait_df),
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
			"grs_or" = numeric(), 
			"grs_lci_or" = numeric(), 
			"grs_uci_or" = numeric(), 
			"n_cases" = integer(),
			"n_controls" = integer(),
			"n_excluded" = integer(),
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

							if (extra_adjustment == "smoking" & trait == "smoker_cases") {
								next()
							}
							#Decide if to run based on the sig analysis since only sensitivity
							decision_sex_group <- gsub("(.)+.eur\\.", "", snp_group)
							decision_sex_group <- gsub("\\.(.)+", "", decision_sex_group)
							decision_snp_group <- paste(grs_trait, ".eur.", decision_sex_group, ".pulit.sig", sep = "")

							if (grs_trait == "whradjbmi") {
								grs_trait <- "res_whr_inv"
							}

							#Make RINTed traits - if WHRadjBMI, use that one directly
                                                        trait_unit_df <- grs_unit_df
                                                        if (grs_trait != "res_whr_inv") {
                                                                rint.formula.exp <- as.formula(paste(grs_trait,
                                                                "~age_assessment+age_squared+dummy_sex+assessment_centre", sep = " "))
                                                                rint_fit_exp <- lm(rint.formula.exp, data = trait_unit_df, na.action = na.exclude)
                                                                trait_unit_df$res_rint_exp <- residuals(rint_fit_exp)
                                                                trait_unit_df$exposure_sd <- qnorm((rank(trait_unit_df$res_rint_exp,
                                                                                na.last = "keep") - 0.5) / sum(!is.na(trait_unit_df$res_rint_exp)))
                                                        } else if (grs_trait == "res_whr_inv") {
                                                                trait_unit_df$exposure_sd <- trait_unit_df$res_whr_inv
                                                        }
							print(warnings())
							n_excluded <- nrow(trait_unit_df[trait_unit_df[, trait] == 2, ])
                	                                trait_df <- trait_unit_df[trait_unit_df[, trait] != 2, ]
                        	                        trait_df[, trait] <- factor(trait_df[, trait])
	
							#Decide which traits to run based on if any sex significant in all.white or brit.irish in the GRS analysis
							if (nrow(trait_df) < 1 | (nrow(grs_results[grs_results$snp_group == decision_snp_group &
                                                                        grs_results$eth_group %in% c("all.white", "brit.irish") &
                                                                        grs_results$sex_group == sex_group &
                                                                        grs_results$grs_unit == grs_unit &
                                                                        grs_results$case_column == trait, ]) < 1)) {
                	                                        next
                        	                        }
							
							#Run regression using Wald method in SD-units and outcome in clinical units
                                	                fit_raw <- run.regression.wald(trait_df, trait, exposure_name = "exposure_sd", outcome_name = trait, extra_adjustment = extra_adjustment)
                                        	        results_df <- current_row_wald(results_df, snp_group, grs_unit, eth_group, sex_group,
                                                	                grs_trait, trait, file, fit_raw, exposure_unit = "sd", outcome_unit = "clin", 
									trait_df, "wald", n_excluded, extra_adjustment)

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

#Clean the phenotype file:
pheno <- qc.pregnant(pheno)
pheno <- qc.below.15.bmi(pheno)

#Read in the results from the logistic regressions
grs_results <- read.table("../results.logistic.regressions.180514/log.results.table.180627.txt", stringsAsFactors = F, 
		header = T, sep = "\t")

comb_groups <- c("bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", "whradjbmi.eur.comb.pulit.sig")
male_groups <- c("bmi.eur.men.pulit.sig", "whr.eur.men.pulit.sig", "whradjbmi.eur.men.pulit.sig")
female_groups <- c("bmi.eur.women.pulit.sig", "whr.eur.women.pulit.sig", "whradjbmi.eur.women.pulit.sig")

adjusted_smoking <- unique(grs_results$case_column[grep("_smoking", grs_results$case_column)])
grs_results <- grs_results[!(grs_results$case_column %in% adjusted_smoking), ]
grs_results <- grs_results[((grs_results$snp_group %in% c(comb_groups) &
                grs_results$sex_group == "comb") |
                (grs_results$snp_group %in% c(male_groups) &
                grs_results$sex_group == "men") |
                (grs_results$snp_group %in% c(female_groups) &
                grs_results$sex_group == "women")) &
                grs_results$grs_unit == "raw_scoresum" &
                grs_results$eth_group %in% c("all.white", "brit.irish"), ]

critical_p <- 0.05/(length(unique(grs_results[!(grs_results$case_column %in% c("t1d_cases_prob", "t2d_cases_prob", "smoker_cases")), 
		"case_column"]))*length(unique(gsub("\\..*", "", grs_results$snp_group))))

for (i in 1:length(comb_groups)) {
	for (grs_unit in unique(grs_results$grs_unit)) {
		for (case_column in unique(grs_results$case_column[!(grs_results$case_column %in% c("t1d_cases_prob", "smoker_cases", "t2d_cases_prob"))])) {
			 if (nrow(grs_results[(grs_results$snp_group %in% comb_groups[i]) &
                                  grs_results$grs_unit == grs_unit & grs_results$case_column == case_column & grs_results$grs_p < critical_p, ]) == 0) {

				grs_results <- grs_results[!(grs_results$snp_group %in% c(comb_groups[i]) &
                                                grs_results$grs_unit == grs_unit & grs_results$case_column == case_column), ]
                        }
                }
        }
}

#Removes the non-significant sex-specific GRSs, except the sensitivity analyses
for (i in 1:length(female_groups)) {
        for (grs_unit in unique(grs_results$grs_unit)) {
                 for (case_column in unique(grs_results$case_column[!(grs_results$case_column %in% c("t1d_cases_prob", "smoker_cases", "t2d_cases_prob"))])) {
                         if (nrow(grs_results[grs_results$snp_group %in% c(female_groups[i], male_groups[i]) &
                                  grs_results$grs_unit == grs_unit & grs_results$case_column == case_column & grs_results$grs_p < critical_p, ]) == 0) {

                                  grs_results <- grs_results[!(grs_results$snp_group %in% c(female_groups[i], male_groups[i]) &
        	                        	grs_results$grs_unit == grs_unit & grs_results$case_column == case_column), ]
                                
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
traits <- c("cad_cases", "stroke_cases", "copd_cases",
                "dementia_cases", "lungcancer_cases", "colorectal_cancer_cases",
                "renal_failure_cases", "breast_cancer_cases",
                "t2d_cases_prob", "t2d_cases_probposs", "t1d_cases_prob",
                "t1d_cases_probposs", 
                "any_infertility_cases",
                "nafld_cases", "cld_cases", "smoker_cases", 
		"aki_cases", "ckd_cases",
		"haem_stroke_cases", "isch_stroke_cases")

results_df <- tsls.cont.mr.grs(pheno, grs_results, snp_groups, eth_groups, grs_unit_groups, sex_groups, traits)

write.table(results_df, "../results.mr.180730/ipd.mr.binary.results.180815.EXTRA.txt", sep = "\t", row.names = F, quote = F)

sink("../results.mr.180730/ipd.mr.binary.results.180815.WARNINGS.txt", append = F)
print(warnings())
sink()

#Calculate heterogeneity values
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
					for (function_name in unique(results_df$function_name)) {

					general <- results_df$eth_group == eth_group &
                                                results_df$trait == trait &
                                                results_df$exposure_unit == exposure_unit &
                                                results_df$outcome_unit == outcome_unit &
						results_df$function_name == function_name

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
}

results_df$cochrans_i2 <- (results_df$cochrans_q - 1)/results_df$cochrans_q
results_df$cochrans_i2[results_df$cochrans_i2 < 0] <- 0

write.table(results_df, "../results.mr.180730/ipd.mr.binary.results.180815.txt", sep = "\t", row.names = F, quote = F)

