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

tsls.cont.mr.grs <- function(pheno, snp_groups, eth_groups, grs_unit_groups, sex_groups, traits) {

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
						for (extra_adjustment in c("-")) {
							grs_trait <- gsub("\\.(.)+", "", snp_group)

							if (extra_adjustment == "smoking" & trait == "smoker_cases") {
								next()
							}
							if (grs_trait == "whradjbmi") {
								grs_trait <- "res_whr_inv"
							}

							grs_sex <- gsub("(.)+\\.eur\\.", "", snp_group)
							grs_sex <- gsub("\\.(.)*", "", grs_sex)
							
							if ((grs_sex != "comb" & grs_sex != sex_group) |
								(trait == "smoker_cases" & length(snp_group[grep("winner", snp_group)]) == 0 ) |
								length(snp_group[grep("winner", snp_group)]) > 0 & 
								!((trait == "t2d_cases_probposs" & grs_trait == "bmi") |
								(trait == "copd_cases" & grs_trait == "whr") |
								(trait == "renal_failure_cases" & grs_trait == "whr") |
								(trait == "smoker_cases" & grs_trait %in% c("bmi", "whr")) |
								(trait == "ckd_cases" & grs_trait %in% c("whr", "res_whr_inv")))){
								next()
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

eth_groups <- c("all.white")

#Get same results for raw scoresum and sd_scoresum
grs_unit_groups <- c("raw_scoresum") 
sex_groups <- c("comb", "women", "men")
traits <- c("t2d_cases_probposs")

types_internal <- c("normal", "half", "randomnoise", "systematicnoise", "unweighted")
types_fdr <- c("0.01.fdr", "0.05.fdr", "0.1.fdr")
x <- paste(rep(c("bmi", "whr", "whradjbmi"), each = length(sex_groups)), "eur", sex_groups, sep = ".")
internal <- paste(rep(x, each = length(types_internal)), "internal", types_internal, sep = ".")
y <- paste(rep(c("bmi", "whr", "whradjbmi"), each = length(c("men", "women"))), "eur", c("men", "women"), sep = ".")
fdr <- paste(rep(y, each = length(types_fdr)), types_fdr, sep = ".")
snp_groups <- c(internal, fdr)

results_df <- tsls.cont.mr.grs(pheno, snp_groups, eth_groups, grs_unit_groups, sex_groups, traits)

write.table(results_df, "../results.mr.180730/sens.ipd.mr.binary.results.180815.EXTRA.txt", sep = "\t", row.names = F, quote = F)

#Comparison of the different weightings for the significant results
traits <- c("t2d_cases_probposs", "cad_cases", "cld_cases", "copd_cases", 
	"lungcancer_cases", "nafld_cases", "renal_failure_cases", 
	"aki_cases", "ckd_cases", "stroke_cases", "isch_stroke_cases", 
	"t1d_cases_probposs", "smoker_cases")
snp_groups <- c(fdr,
		"bmi.eur.comb.pulit.sig", "whr.eur.comb.pulit.sig", 
		"whradjbmi.eur.comb.pulit.sig",
		"bmi.eur.men.pulit.sig", "bmi.eur.women.pulit.sig",
                "whr.eur.men.pulit.sig", "whr.eur.women.pulit.sig",
                "whradjbmi.eur.men.pulit.sig", "whradjbmi.eur.women.pulit.sig", 
		"bmi.eur.men.pulit.phet", "bmi.eur.women.pulit.phet",
                "whr.eur.men.pulit.phet", "whr.eur.women.pulit.phet",
                "whradjbmi.eur.men.pulit.phet", "whradjbmi.eur.women.pulit.phet",
		"bmi.eur.men.internal.unweighted", "bmi.eur.women.internal.unweighted",
                "whr.eur.men.internal.unweighted", "whr.eur.women.internal.unweighted",
                "whradjbmi.eur.men.internal.unweighted", "whradjbmi.eur.women.internal.unweighted", 
		"bmi.eur.comb.internal.unweighted", "whr.eur.comb.internal.unweighted", 
		"whradjbmi.eur.comb.internal.unweighted", 
		"bmi.eur.men.pulit", "whr.eur.men.pulit", "whradjbmi.eur.men.pulit", 
		"bmi.eur.women.pulit", "whr.eur.women.pulit", "whradjbmi.eur.women.pulit", 
		"bmi.eur.comb.pulit.winner", "whr.eur.comb.pulit.winner",
                "whradjbmi.eur.comb.pulit.winner",
                "bmi.eur.men.pulit.winner", "whr.eur.men.pulit.winner",
                "whradjbmi.eur.men.pulit.winner",
                "bmi.eur.women.pulit.winner", "whr.eur.women.pulit.winner",
                "whradjbmi.eur.women.pulit.winner", 
		"bmi.eur.comb.pulit.winner_unweighted", "whr.eur.comb.pulit.winner_unweighted", 
		"whradjbmi.eur.comb.pulit.winner_unweighted")

results_df_all_traits <- tsls.cont.mr.grs(pheno, snp_groups, eth_groups, grs_unit_groups, sex_groups, traits)

results_df <- rbind(results_df, results_df_all_traits)
results_df <- results_df[!duplicated(results_df[]), ]

sink("../results.mr.180730/sens.ipd.mr.binary.results.180815.WARNINGS.txt", append = F)
print(warnings())
sink()

#Calculate heterogeneity values
results_df$cochran_weights <- 1 / (results_df$grs_se^2)

snp_groups <- c(snp_groups, internal)
snp_groups <- snp_groups[!duplicated(snp_groups)]
comb_snp_groups <- snp_groups[grep("\\.comb\\.", snp_groups)]
women_snp_groups <- snp_groups[grep("\\.women\\.", snp_groups)]
men_snp_groups <- snp_groups[grep("\\.men\\.", snp_groups)]

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

write.table(results_df, "../results.mr.180730/sens.ipd.mr.binary.results.180815.txt", sep = "\t", row.names = F, quote = F)

